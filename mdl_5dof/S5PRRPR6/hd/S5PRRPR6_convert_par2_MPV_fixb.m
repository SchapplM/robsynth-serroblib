% Return the minimum parameter vector for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [22x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5PRRPR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t123 = 2 * pkin(8) * mrSges(6,3) + Ifges(6,2);
t115 = sin(pkin(10));
t116 = cos(pkin(10));
t122 = t115 * t116;
t121 = Ifges(5,4) * t122;
t120 = -pkin(8) * m(6) - mrSges(6,3);
t118 = (pkin(4) ^ 2);
t119 = t118 * m(6) + Ifges(4,2) + Ifges(5,3);
t117 = pkin(8) ^ 2;
t113 = t116 ^ 2;
t112 = t115 ^ 2;
t111 = t120 * pkin(4) + Ifges(5,5);
t110 = m(6) * t117 + Ifges(5,1) + t123;
t109 = Ifges(5,2) + (t117 + t118) * m(6) + t123;
t1 = [m(2) + m(3) + m(4); Ifges(3,3) + 2 * pkin(7) * mrSges(4,3) + (pkin(2) ^ 2 + pkin(7) ^ 2) * m(4) + t119; m(4) * pkin(2) + mrSges(3,1); -pkin(7) * m(4) + mrSges(3,2) - mrSges(4,3); t112 * t109 + t113 * t110 + Ifges(4,1) - t119 - 0.2e1 * t121; t115 * Ifges(5,6) - t116 * t111 + Ifges(4,4); Ifges(4,5) + (t113 - t112) * Ifges(5,4) + (-t109 + t110) * t122; -t116 * Ifges(5,6) - t115 * t111 + Ifges(4,6); t113 * t109 + t112 * t110 + Ifges(4,3) + 0.2e1 * t121; mrSges(4,1); mrSges(4,2); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t120; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
