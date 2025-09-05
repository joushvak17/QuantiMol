import sys
from PySide6.QtWidgets import QApplication


def setup_application() -> QApplication:
    """Setup the QApplication with necessary settings."""
    app = QApplication(sys.argv)

    # Set application properties
    app.setApplicationName("QuantiMol")
    app.setApplicationVersion("0.0.1")
    app.setOrganizationName("joushvak")
    app.setOrganizationDomain("https://github.com/joushvak17")

    # TODO: Set application icon
    # app_icon = QIcon("gui/resources/images/Icon.ico")
    # app.setWindowIcon(app_icon)

    # TODO: Set application style (optional)
    # app.setStyle("")

    # TODO: Look into enabling high DPI support
    # app.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    # app.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    return app


def main() -> int:
    try:
        app = setup_application()
        
        # TODO: This should implement a QMainWindow

        exit_code = app.exec()

        return exit_code

    except Exception as e:
        print(f"Failed to start the application: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
