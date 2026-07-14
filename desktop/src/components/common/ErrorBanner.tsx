interface Props {
  message: string;
}

export function ErrorBanner({ message }: Props) {
  return (
    <div className="error-banner">
      <strong>Error:</strong> {message}
    </div>
  );
}
