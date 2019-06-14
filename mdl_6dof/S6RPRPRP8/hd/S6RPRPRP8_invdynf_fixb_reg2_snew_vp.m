% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RPRPRP8_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:03:28
% EndTime: 2019-05-05 18:03:33
% DurationCPUTime: 5.90s
% Computational Cost: add. (14123->277), mult. (31245->295), div. (0->0), fcn. (21193->8), ass. (0->173)
t1932 = cos(qJ(3));
t1975 = qJD(1) * t1932;
t1956 = qJD(3) * t1975;
t1929 = sin(qJ(3));
t1961 = t1929 * qJDD(1);
t1904 = -t1956 - t1961;
t1976 = qJD(1) * t1929;
t1958 = qJD(3) * t1976;
t1960 = t1932 * qJDD(1);
t1905 = -t1958 + t1960;
t1925 = sin(pkin(9));
t1926 = cos(pkin(9));
t1952 = -t1926 * t1904 + t1925 * t1905;
t1950 = qJDD(5) + t1952;
t1900 = -t1925 * t1976 + t1926 * t1975;
t1928 = sin(qJ(5));
t1931 = cos(qJ(5));
t1881 = -t1931 * qJD(3) + t1928 * t1900;
t1883 = t1928 * qJD(3) + t1931 * t1900;
t1970 = t1883 * t1881;
t1836 = t1950 + t1970;
t1880 = t1883 ^ 2;
t1898 = (t1925 * t1932 + t1926 * t1929) * qJD(1);
t1896 = qJD(5) + t1898;
t1981 = t1896 ^ 2;
t1986 = -t1880 - t1981;
t1817 = t1931 * t1836 + t1928 * t1986;
t1876 = t1925 * t1904 + t1926 * t1905;
t1944 = -t1881 * qJD(5) + t1928 * qJDD(3) + t1931 * t1876;
t1971 = t1881 * t1896;
t1940 = t1944 - t1971;
t1794 = t1925 * t1817 + t1926 * t1940;
t1796 = t1926 * t1817 - t1925 * t1940;
t1775 = t1932 * t1794 + t1929 * t1796;
t1815 = t1928 * t1836 - t1931 * t1986;
t1930 = sin(qJ(1));
t1933 = cos(qJ(1));
t2008 = t1930 * t1775 + t1933 * t1815;
t2007 = t1933 * t1775 - t1930 * t1815;
t1781 = t1929 * t1794 - t1932 * t1796;
t1831 = t1944 + t1971;
t1953 = -t1931 * qJDD(3) + t1928 * t1876;
t1945 = (-qJD(5) + t1896) * t1883 - t1953;
t1984 = -t1931 * t1831 + t1928 * t1945;
t1860 = t1881 ^ 2;
t1842 = t1880 + t1860;
t1983 = t1928 * t1831 + t1931 * t1945;
t1995 = t1926 * t1842 + t1925 * t1983;
t1996 = -t1925 * t1842 + t1926 * t1983;
t2002 = t1929 * t1996 + t1932 * t1995;
t2006 = t1930 * t2002 + t1933 * t1984;
t2005 = t1930 * t1984 - t1933 * t2002;
t2001 = -t1929 * t1995 + t1932 * t1996;
t1837 = t1950 - t1970;
t1985 = -t1981 - t1860;
t1989 = -t1928 * t1837 + t1931 * t1985;
t2000 = t1925 * t1989;
t1999 = t1926 * t1989;
t1990 = t1931 * t1837 + t1928 * t1985;
t1998 = t1930 * t1990;
t1997 = t1933 * t1990;
t1922 = t1929 ^ 2;
t1935 = qJD(1) ^ 2;
t1947 = qJD(3) * pkin(3) - qJ(4) * t1975;
t1914 = -t1933 * g(1) - t1930 * g(2);
t1948 = -qJDD(1) * qJ(2) - t1914;
t1977 = pkin(7) + pkin(1);
t1982 = t1904 * pkin(3) + (t1922 * qJ(4) + t1977) * t1935 - t1947 * t1975 - qJDD(4) + t1948;
t1980 = t1898 ^ 2;
t1979 = t1900 ^ 2;
t1978 = 2 * qJD(4);
t1974 = qJD(2) * qJD(1);
t1973 = qJD(3) * t1898;
t1972 = qJD(3) * t1900;
t1969 = t1900 * t1898;
t1968 = t1922 * t1935;
t1963 = t1932 * t1935;
t1913 = t1930 * g(1) - t1933 * g(2);
t1943 = -t1935 * qJ(2) + qJDD(2) - t1913;
t1939 = -t1977 * qJDD(1) + t1943;
t1878 = -t1932 * g(3) + t1929 * t1939;
t1856 = -pkin(3) * t1968 + t1904 * qJ(4) - qJD(3) * t1947 + t1878;
t1938 = t1932 * t1939;
t1937 = t1938 - t1905 * qJ(4) + qJDD(3) * pkin(3) + (-qJD(3) * qJ(4) * qJD(1) - pkin(3) * t1963 + g(3)) * t1929;
t1824 = t1926 * t1856 - t1898 * t1978 + t1925 * t1937;
t1870 = t1898 * pkin(4) - t1900 * pkin(8);
t1934 = qJD(3) ^ 2;
t1810 = -t1934 * pkin(4) + qJDD(3) * pkin(8) - t1898 * t1870 + t1824;
t1861 = t1952 + t1972;
t1955 = -t1876 + t1973;
t1936 = t1861 * pkin(4) + t1955 * pkin(8) + 0.2e1 * t1974 - t1982;
t1787 = t1931 * t1810 + t1928 * t1936;
t1923 = t1932 ^ 2;
t1962 = t1922 + t1923;
t1959 = -0.2e1 * t1974;
t1957 = t1929 * t1963;
t1786 = -t1928 * t1810 + t1931 * t1936;
t1954 = t1925 * t1856 - t1926 * t1937;
t1946 = -t1883 * qJD(5) - t1953;
t1942 = t1948 + t1959;
t1809 = -qJDD(3) * pkin(4) - t1934 * pkin(8) + (t1978 + t1870) * t1900 + t1954;
t1917 = -t1923 * t1935 - t1934;
t1916 = -t1934 - t1968;
t1912 = -qJDD(3) - t1957;
t1911 = qJDD(3) - t1957;
t1910 = t1962 * t1935;
t1909 = t1930 * qJDD(1) + t1933 * t1935;
t1908 = t1933 * qJDD(1) - t1930 * t1935;
t1907 = t1962 * qJDD(1);
t1906 = -0.2e1 * t1958 + t1960;
t1903 = 0.2e1 * t1956 + t1961;
t1897 = qJDD(1) * pkin(1) - t1943;
t1892 = t1935 * pkin(1) + t1942;
t1890 = -t1934 - t1979;
t1889 = t1977 * t1935 + t1942;
t1887 = t1932 * t1912 - t1929 * t1917;
t1886 = -t1929 * t1911 + t1932 * t1916;
t1885 = t1929 * t1912 + t1932 * t1917;
t1884 = t1932 * t1911 + t1929 * t1916;
t1877 = t1929 * g(3) + t1938;
t1874 = -qJDD(3) - t1969;
t1873 = qJDD(3) - t1969;
t1871 = -t1934 - t1980;
t1864 = -t1876 - t1973;
t1862 = -t1952 + t1972;
t1859 = -t1979 - t1980;
t1858 = t1881 * pkin(5) - t1883 * qJ(6);
t1857 = t1959 + t1982;
t1853 = t1926 * t1874 - t1925 * t1890;
t1852 = t1925 * t1874 + t1926 * t1890;
t1850 = -t1929 * t1877 + t1932 * t1878;
t1849 = t1932 * t1877 + t1929 * t1878;
t1840 = t1926 * t1871 - t1925 * t1873;
t1839 = t1925 * t1871 + t1926 * t1873;
t1835 = t1926 * t1862 - t1925 * t1864;
t1834 = t1925 * t1862 + t1926 * t1864;
t1828 = t1896 * t1883 - t1946;
t1827 = (qJD(5) + t1896) * t1883 + t1953;
t1826 = -t1929 * t1852 + t1932 * t1853;
t1825 = t1932 * t1852 + t1929 * t1853;
t1823 = -0.2e1 * qJD(4) * t1900 - t1954;
t1812 = -t1929 * t1839 + t1932 * t1840;
t1811 = t1932 * t1839 + t1929 * t1840;
t1807 = -t1929 * t1834 + t1932 * t1835;
t1806 = t1932 * t1834 + t1929 * t1835;
t1801 = t1925 * t1828 + t1999;
t1799 = -t1926 * t1828 + t2000;
t1797 = t1925 * t1827 + t1999;
t1795 = -t1926 * t1827 + t2000;
t1793 = -t1925 * t1823 + t1926 * t1824;
t1792 = t1926 * t1823 + t1925 * t1824;
t1785 = -t1946 * pkin(5) + (pkin(5) * t1896 - (2 * qJD(6))) * t1883 + t1809 - t1940 * qJ(6);
t1784 = -t1950 * pkin(5) - qJ(6) * t1981 + t1883 * t1858 + qJDD(6) - t1786;
t1783 = -pkin(5) * t1981 + t1950 * qJ(6) + 0.2e1 * qJD(6) * t1896 - t1881 * t1858 + t1787;
t1782 = -t1929 * t1799 + t1932 * t1801;
t1780 = t1932 * t1799 + t1929 * t1801;
t1778 = -t1929 * t1795 + t1932 * t1797;
t1776 = t1932 * t1795 + t1929 * t1797;
t1774 = -t1929 * t1792 + t1932 * t1793;
t1773 = t1932 * t1792 + t1929 * t1793;
t1768 = -t1928 * t1786 + t1931 * t1787;
t1767 = t1931 * t1786 + t1928 * t1787;
t1766 = t1926 * t1768 + t1925 * t1809;
t1765 = t1925 * t1768 - t1926 * t1809;
t1764 = t1931 * t1783 + t1928 * t1784;
t1763 = t1928 * t1783 - t1931 * t1784;
t1762 = t1926 * t1764 + t1925 * t1785;
t1761 = t1925 * t1764 - t1926 * t1785;
t1760 = -t1929 * t1765 + t1932 * t1766;
t1759 = t1932 * t1765 + t1929 * t1766;
t1758 = -t1929 * t1761 + t1932 * t1762;
t1757 = t1932 * t1761 + t1929 * t1762;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1909, -t1908, 0, -t1930 * t1913 + t1933 * t1914, 0, 0, 0, 0, 0, 0, 0, t1909, t1908, -t1933 * t1892 - t1930 * t1897, 0, 0, 0, 0, 0, 0, t1930 * t1884 + t1933 * t1903, t1930 * t1885 + t1933 * t1906, -t1930 * t1907 - t1933 * t1910, t1930 * t1849 - t1933 * t1889, 0, 0, 0, 0, 0, 0, t1930 * t1811 + t1933 * t1861, t1930 * t1825 - t1933 * t1955, t1930 * t1806 + t1933 * t1859, t1930 * t1773 - t1933 * t1857, 0, 0, 0, 0, 0, 0, t1930 * t1776 + t1997, -t2008, t2006, t1930 * t1759 + t1933 * t1767, 0, 0, 0, 0, 0, 0, t1930 * t1780 + t1997, t2006, t2008, t1930 * t1757 + t1933 * t1763; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1908, -t1909, 0, t1933 * t1913 + t1930 * t1914, 0, 0, 0, 0, 0, 0, 0, -t1908, t1909, -t1930 * t1892 + t1933 * t1897, 0, 0, 0, 0, 0, 0, -t1933 * t1884 + t1930 * t1903, -t1933 * t1885 + t1930 * t1906, t1933 * t1907 - t1930 * t1910, -t1933 * t1849 - t1930 * t1889, 0, 0, 0, 0, 0, 0, -t1933 * t1811 + t1930 * t1861, -t1933 * t1825 - t1930 * t1955, -t1933 * t1806 + t1930 * t1859, -t1933 * t1773 - t1930 * t1857, 0, 0, 0, 0, 0, 0, -t1933 * t1776 + t1998, t2007, t2005, -t1933 * t1759 + t1930 * t1767, 0, 0, 0, 0, 0, 0, -t1933 * t1780 + t1998, t2005, -t2007, -t1933 * t1757 + t1930 * t1763; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1886, t1887, 0, t1850, 0, 0, 0, 0, 0, 0, t1812, t1826, t1807, t1774, 0, 0, 0, 0, 0, 0, t1778, t1781, t2001, t1760, 0, 0, 0, 0, 0, 0, t1782, t2001, -t1781, t1758; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1935, -qJDD(1), 0, t1914, 0, 0, 0, 0, 0, 0, 0, t1935, qJDD(1), -t1892, 0, 0, 0, 0, 0, 0, t1903, t1906, -t1910, -t1889, 0, 0, 0, 0, 0, 0, t1861, -t1955, t1859, -t1857, 0, 0, 0, 0, 0, 0, t1990, -t1815, t1984, t1767, 0, 0, 0, 0, 0, 0, t1990, t1984, t1815, t1763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1935, 0, t1913, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t1935, t1897, 0, 0, 0, 0, 0, 0, -t1884, -t1885, t1907, -t1849, 0, 0, 0, 0, 0, 0, -t1811, -t1825, -t1806, -t1773, 0, 0, 0, 0, 0, 0, -t1776, t1775, -t2002, -t1759, 0, 0, 0, 0, 0, 0, -t1780, -t2002, -t1775, -t1757; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1886, t1887, 0, t1850, 0, 0, 0, 0, 0, 0, t1812, t1826, t1807, t1774, 0, 0, 0, 0, 0, 0, t1778, t1781, t2001, t1760, 0, 0, 0, 0, 0, 0, t1782, t2001, -t1781, t1758; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1886, t1887, 0, t1850, 0, 0, 0, 0, 0, 0, t1812, t1826, t1807, t1774, 0, 0, 0, 0, 0, 0, t1778, t1781, t2001, t1760, 0, 0, 0, 0, 0, 0, t1782, t2001, -t1781, t1758; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1935, -qJDD(1), t1892, 0, 0, 0, 0, 0, 0, -t1903, -t1906, t1910, t1889, 0, 0, 0, 0, 0, 0, -t1861, t1955, -t1859, t1857, 0, 0, 0, 0, 0, 0, -t1990, t1815, -t1984, -t1767, 0, 0, 0, 0, 0, 0, -t1990, -t1984, -t1815, -t1763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1935, -t1897, 0, 0, 0, 0, 0, 0, t1884, t1885, -t1907, t1849, 0, 0, 0, 0, 0, 0, t1811, t1825, t1806, t1773, 0, 0, 0, 0, 0, 0, t1776, -t1775, t2002, t1759, 0, 0, 0, 0, 0, 0, t1780, t2002, t1775, t1757; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1916, t1912, -t1961, t1878, 0, 0, 0, 0, 0, 0, t1840, t1853, t1835, t1793, 0, 0, 0, 0, 0, 0, t1797, -t1796, t1996, t1766, 0, 0, 0, 0, 0, 0, t1801, t1996, t1796, t1762; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1911, t1917, -t1960, t1877, 0, 0, 0, 0, 0, 0, t1839, t1852, t1834, t1792, 0, 0, 0, 0, 0, 0, t1795, -t1794, t1995, t1765, 0, 0, 0, 0, 0, 0, t1799, t1995, t1794, t1761; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1903, t1906, -t1910, -t1889, 0, 0, 0, 0, 0, 0, t1861, -t1955, t1859, -t1857, 0, 0, 0, 0, 0, 0, t1990, -t1815, t1984, t1767, 0, 0, 0, 0, 0, 0, t1990, t1984, t1815, t1763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1871, t1874, t1862, t1824, 0, 0, 0, 0, 0, 0, t1989, -t1817, t1983, t1768, 0, 0, 0, 0, 0, 0, t1989, t1983, t1817, t1764; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1873, t1890, t1864, t1823, 0, 0, 0, 0, 0, 0, -t1827, -t1940, t1842, -t1809, 0, 0, 0, 0, 0, 0, -t1828, t1842, t1940, -t1785; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1861, -t1955, t1859, -t1857, 0, 0, 0, 0, 0, 0, t1990, -t1815, t1984, t1767, 0, 0, 0, 0, 0, 0, t1990, t1984, t1815, t1763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1985, -t1836, t1945, t1787, 0, 0, 0, 0, 0, 0, t1985, t1945, t1836, t1783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1837, t1986, -t1831, t1786, 0, 0, 0, 0, 0, 0, t1837, -t1831, -t1986, -t1784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1827, t1940, -t1842, t1809, 0, 0, 0, 0, 0, 0, t1828, -t1842, -t1940, t1785; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1985, t1945, t1836, t1783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1828, -t1842, -t1940, t1785; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1837, t1831, t1986, t1784;];
f_new_reg  = t1;