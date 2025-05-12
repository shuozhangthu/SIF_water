import numpy as np
import scipy.io as scio
from thundersvm import SVC  # GPU accelerated


def load_climate_data(data_path):
    """Load climate data with precipitation and temperature features
    Args:
        data_path: path to .mat data file
    Returns:
        Tuple of (features, labels) where:
        - features: ndarray of shape (n_samples, 2) [precipitation, temperature]
        - labels: ndarray of shape (n_samples,) [DOY_SIF_PAR]
    """
    mat_data = scio.loadmat(data_path)
    climate_data = mat_data['var']

    # First two columns are features [precipitation, temperature]
    features = climate_data[:, :2]
    # Last column is label [DOY_SIF_PAR]
    labels = climate_data[:, -1]

    return features, labels


def split_data(features, labels, test_size=0.2):

    n_samples = features.shape[0]
    split_idx = int(n_samples * (1 - test_size))

    X_train, X_test = features[:split_idx], features[split_idx:]
    y_train, y_test = labels[:split_idx], labels[split_idx:]

    return X_train, X_test, y_train, y_test


def train_and_evaluate(X_train, X_test, y_train, y_test):
    """Train SVM model and evaluate performance"""
    # Initialize GPU-accelerated SVM
    model = SVC(kernel='linear', C=1.0, gamma='auto')
    model.fit(X_train, y_train)

    print("\nModel Coefficients:")
    print(f"Precipitation coef: {model.coef_[0][0]:.4f}")
    print(f"Temperature coef: {model.coef_[0][1]:.4f}")
    print(f"Intercept: {model.intercept_[0]:.4f}")

    accuracy = model.score(X_test, y_test)
    print(f"\nTest Accuracy: {accuracy:.4f}")

    return model


def main():
    DATA_PATH = "./climate_data.mat"

    print("Loading climate data...")
    features, labels = load_climate_data(DATA_PATH)

    print("Splitting dataset...")
    X_train, X_test, y_train, y_test = split_data(features, labels)

    print("Training SVM model...")
    model = train_and_evaluate(X_train, X_test, y_train, y_test)


if __name__ == "__main__":
    main()


