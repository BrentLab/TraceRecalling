%ifarch x86_64
%define estbin est_genome_ambig.x86_64
%else
%define estbin est_genome_ambig.i386
%endif

Summary: Trace Recalling
Name: trace_recalling
Version: 0.5
Release: 2
License: ?
Group: ?
URL: ?
Source0: %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
Requires: eval-common

%description

%prep
%setup -q

%build
cp %{estbin} est_genome_ambig
%{__perl} Makefile.PL INSTALLDIRS="vendor" PREFIX="%{buildroot}%{_prefix}"


%install
rm -rf $RPM_BUILD_ROOT
make install

find $RPM_BUILD_ROOT \( -name perllocal.pod -o -name .packlist \) -exec rm -v {} \;

rm -rf $RPM_BUILD_ROOT/usr/share/man

find $RPM_BUILD_ROOT/usr -type f -print | \
    sed "s@^$RPM_BUILD_ROOT@@g" | \
    grep -v perllocal.pod | \
    grep -v "\.packlist" > trace_recalling-%{version}-filelist




%clean
rm -rf $RPM_BUILD_ROOT

%files -f trace_recalling-%{version}-filelist
%defattr(-,root,root,-)
%doc README

%changelog
* Wed Feb 28 2007 Brian Koebbe <koebbe@ardor.mblab.wustl.edu> - 0.5-2
- a few minor path changes

* Fri Feb  9 2007 Aaron Tenney <tenney@ardor.mblab.wustl.edu> - 0.5-1
- Initial build.

