% Return the minimum parameter vector for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% MPV [27x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPRR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t120 = (m(6) + m(7));
t123 = (pkin(8) ^ 2);
t125 = (pkin(4) ^ 2);
t124 = (pkin(5) ^ 2);
t135 = (t124 * m(7) + Ifges(6,2));
t130 = 2 * pkin(8) * mrSges(6,3) + t135;
t106 = Ifges(5,2) + (t123 + t125) * t120 + t130;
t107 = t123 * t120 + Ifges(5,1) + t130;
t136 = -t106 + t107;
t134 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t116 = sin(pkin(11));
t118 = cos(pkin(11));
t133 = t116 * t118;
t112 = t116 ^ 2;
t113 = t118 ^ 2;
t132 = t113 - t112;
t131 = Ifges(5,4) * t133;
t129 = pkin(9) * m(7) + mrSges(7,3);
t128 = -pkin(8) * t120 - mrSges(6,3);
t109 = pkin(4) * t120 + mrSges(5,1);
t127 = -t116 * mrSges(5,2) + t118 * t109;
t110 = -pkin(7) * m(4) + mrSges(3,2) - mrSges(4,3);
t111 = m(4) * pkin(2) + mrSges(3,1);
t117 = sin(pkin(10));
t119 = cos(pkin(10));
t126 = -t117 * t110 + t119 * t111;
t122 = pkin(9) ^ 2;
t108 = pkin(4) * t128 + Ifges(5,5);
t1 = [Ifges(2,3) + Ifges(3,3) + Ifges(4,2) + t112 * t107 + 0.2e1 * t131 + t113 * t106 + (2 * pkin(7) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * m(4)) + 0.2e1 * t126 * pkin(1); mrSges(2,1) + t126; t119 * t110 + t117 * t111 + mrSges(2,2); m(3) + m(4); t136 * t132 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t131; t132 * Ifges(5,4) + t136 * t133 + Ifges(4,4); -t116 * Ifges(5,6) + t118 * t108 + Ifges(4,5); t118 * Ifges(5,6) + t116 * t108 + Ifges(4,6); 0.2e1 * pkin(3) * t127 + (t125 * t120) + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t127; t118 * mrSges(5,2) + t116 * t109 + mrSges(4,2); mrSges(5,3) - t128; m(5) + t120; m(7) * t122 + Ifges(6,1) + t134 - t135; pkin(5) * t129 + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t122 + t124) * m(7) + t134; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t129; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
