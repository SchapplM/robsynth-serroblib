% Return the minimum parameter vector for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t121 = (pkin(9) ^ 2);
t133 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t106 = m(7) * t121 + Ifges(6,1) + t133;
t123 = (pkin(5) ^ 2);
t110 = t123 * m(7) + Ifges(6,2);
t134 = t106 - t110;
t120 = (m(4) + m(5));
t116 = sin(pkin(11));
t118 = cos(pkin(11));
t132 = t116 * t118;
t113 = t116 ^ 2;
t114 = t118 ^ 2;
t131 = t114 - t113;
t130 = -pkin(8) * m(5) - mrSges(5,3);
t129 = pkin(9) * m(7) + mrSges(7,3);
t107 = t129 * pkin(5) + Ifges(6,4);
t128 = t107 * t132;
t127 = (mrSges(4,3) - t130);
t105 = -pkin(7) * t120 + mrSges(3,2) - t127;
t108 = pkin(2) * t120 + mrSges(3,1);
t117 = sin(pkin(10));
t119 = cos(pkin(10));
t126 = -t117 * t105 + t119 * t108;
t109 = mrSges(6,2) - t129;
t111 = m(7) * pkin(5) + mrSges(6,1);
t125 = -t116 * t109 + t118 * t111;
t124 = pkin(3) ^ 2;
t122 = pkin(8) ^ 2;
t112 = t122 + t124;
t1 = [Ifges(2,3) + Ifges(3,3) + Ifges(4,2) + Ifges(5,2) + t113 * t106 + 0.2e1 * t128 + t114 * t110 + (2 * pkin(8) * mrSges(5,3)) + t112 * m(5) + (2 * pkin(7) * t127) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * t120) + 0.2e1 * t126 * pkin(1); mrSges(2,1) + t126; t119 * t105 + t117 * t108 + mrSges(2,2); m(3) + t120; Ifges(4,1) - Ifges(4,2) + (-t112 + t122) * m(5); Ifges(4,4); t130 * pkin(3) + Ifges(4,5); Ifges(4,6); t124 * m(5) + Ifges(4,3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); t134 * t131 + Ifges(5,1) - Ifges(5,2) - 0.4e1 * t128; t131 * t107 + t134 * t132 + Ifges(5,4); t118 * Ifges(6,5) - t116 * Ifges(6,6) + Ifges(5,5); t116 * Ifges(6,5) + t118 * Ifges(6,6) + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t121 + t123) * m(7)) + 0.2e1 * t125 * pkin(4) + t133; mrSges(5,1) + t125; t118 * t109 + t116 * t111 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
