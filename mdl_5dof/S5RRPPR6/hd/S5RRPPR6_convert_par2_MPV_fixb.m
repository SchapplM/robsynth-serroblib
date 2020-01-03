% Return the minimum parameter vector for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% MPV [23x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPPR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t116 = (pkin(7) ^ 2);
t117 = (pkin(4) ^ 2);
t125 = 2 * pkin(7) * mrSges(6,3) + Ifges(6,2);
t103 = Ifges(5,2) + (t116 + t117) * m(6) + t125;
t104 = m(6) * t116 + Ifges(5,1) + t125;
t112 = sin(pkin(9));
t107 = t112 ^ 2;
t114 = cos(pkin(9));
t109 = t114 ^ 2;
t124 = t112 * t114;
t121 = Ifges(5,4) * t124;
t100 = t107 * t103 + t109 * t104 + Ifges(4,1) - 0.2e1 * t121;
t106 = t117 * m(6) + Ifges(4,2) + Ifges(5,3);
t126 = t100 - t106;
t113 = sin(pkin(8));
t115 = cos(pkin(8));
t123 = t113 * t115;
t108 = t113 ^ 2;
t110 = t115 ^ 2;
t122 = t110 - t108;
t120 = -pkin(7) * m(6) - mrSges(6,3);
t105 = t120 * pkin(4) + Ifges(5,5);
t102 = t112 * Ifges(5,6) - t114 * t105 + Ifges(4,4);
t119 = t102 * t123;
t118 = t115 * mrSges(4,1) - t113 * mrSges(4,2);
t101 = -t114 * Ifges(5,6) - t112 * t105 + Ifges(4,6);
t99 = Ifges(4,5) + (t109 - t107) * Ifges(5,4) + (-t103 + t104) * t124;
t1 = [Ifges(2,3) + Ifges(3,2) + t108 * t100 + 0.2e1 * t119 + t110 * t106 + (2 * pkin(6) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(6) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(6) * m(3) + mrSges(2,2) - mrSges(3,3); t126 * t122 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t119; t122 * t102 + t126 * t123 + Ifges(3,4); -t113 * t101 + t115 * t99 + Ifges(3,5); t115 * t101 + t113 * t99 + Ifges(3,6); 0.2e1 * pkin(2) * t118 + t109 * t103 + t107 * t104 + Ifges(3,3) + Ifges(4,3) + 0.2e1 * t121; mrSges(3,1) + t118; t113 * mrSges(4,1) + t115 * mrSges(4,2) + mrSges(3,2); mrSges(4,3); m(4); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t120; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
