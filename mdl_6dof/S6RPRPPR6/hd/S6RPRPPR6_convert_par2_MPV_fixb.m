% Return the minimum parameter vector for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [26x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPPR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t122 = (pkin(8) ^ 2);
t123 = (pkin(5) ^ 2);
t131 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t109 = Ifges(6,2) + (t122 + t123) * m(7) + t131;
t110 = m(7) * t122 + Ifges(6,1) + t131;
t118 = sin(pkin(10));
t113 = t118 ^ 2;
t120 = cos(pkin(10));
t115 = t120 ^ 2;
t130 = t118 * t120;
t127 = Ifges(6,4) * t130;
t106 = t109 * t113 + t110 * t115 + Ifges(5,1) - 0.2e1 * t127;
t112 = m(7) * t123 + Ifges(5,2) + Ifges(6,3);
t133 = t106 - t112;
t132 = (pkin(7) * m(4));
t119 = sin(pkin(9));
t121 = cos(pkin(9));
t129 = t119 * t121;
t114 = t119 ^ 2;
t116 = t121 ^ 2;
t128 = t116 - t114;
t126 = -m(7) * pkin(8) - mrSges(7,3);
t111 = pkin(5) * t126 + Ifges(6,5);
t108 = Ifges(6,6) * t118 - t111 * t120 + Ifges(5,4);
t125 = t108 * t129;
t124 = mrSges(5,1) * t121 - mrSges(5,2) * t119;
t107 = -Ifges(6,6) * t120 - t111 * t118 + Ifges(5,6);
t105 = Ifges(5,5) + (t115 - t113) * Ifges(6,4) + (-t109 + t110) * t130;
t1 = [0.2e1 * t125 + t114 * t106 + t116 * t112 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + ((2 * mrSges(4,3) + t132) * pkin(7)); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t132; mrSges(3,3); m(3) + m(4); t128 * t133 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t125; t108 * t128 + t129 * t133 + Ifges(4,4); t105 * t121 - t107 * t119 + Ifges(4,5); t105 * t119 + t107 * t121 + Ifges(4,6); 0.2e1 * pkin(3) * t124 + t115 * t109 + t113 * t110 + Ifges(4,3) + Ifges(5,3) + 0.2e1 * t127; mrSges(4,1) + t124; mrSges(5,1) * t119 + mrSges(5,2) * t121 + mrSges(4,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t126; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
