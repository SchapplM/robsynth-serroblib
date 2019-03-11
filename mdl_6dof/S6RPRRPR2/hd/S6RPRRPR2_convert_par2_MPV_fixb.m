% Return the minimum parameter vector for
% S6RPRRPR2
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t127 = (m(4) + m(5));
t140 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t123 = sin(pkin(11));
t125 = cos(pkin(11));
t139 = t123 * t125;
t138 = Ifges(6,4) * t139;
t128 = pkin(9) ^ 2;
t130 = (pkin(5) ^ 2);
t112 = Ifges(6,2) + (t128 + t130) * m(7) + t140;
t114 = m(7) * t128 + Ifges(6,1) + t140;
t119 = t123 ^ 2;
t120 = t125 ^ 2;
t137 = t120 * t112 + t119 * t114 + Ifges(5,2);
t136 = pkin(8) * m(5) + mrSges(5,3);
t135 = -pkin(9) * m(7) - mrSges(7,3);
t134 = (2 * pkin(8) * mrSges(5,3)) + t137 + 0.2e1 * t138;
t118 = m(7) * pkin(5) + mrSges(6,1);
t133 = -t123 * mrSges(6,2) + t125 * t118;
t116 = -pkin(7) * t127 + mrSges(3,2) - mrSges(4,3);
t117 = pkin(2) * t127 + mrSges(3,1);
t124 = sin(pkin(10));
t126 = cos(pkin(10));
t132 = -t124 * t116 + t126 * t117;
t131 = pkin(3) ^ 2;
t129 = pkin(8) ^ 2;
t115 = pkin(5) * t135 + Ifges(6,5);
t1 = [Ifges(2,3) + Ifges(3,3) + t131 * m(5) + Ifges(4,2) + (2 * pkin(7) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * t127) + 0.2e1 * t132 * pkin(1); mrSges(2,1) + t132; t126 * t116 + t124 * t117 + mrSges(2,2); m(3) + t127; Ifges(4,1) - Ifges(4,2) + (t129 - t131) * m(5) + t134; t136 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t129 + t131) * m(5) + t134; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t136; t119 * t112 + t120 * t114 + Ifges(5,1) - t137 - 0.4e1 * t138; Ifges(5,4) + (t120 - t119) * Ifges(6,4) + (-t112 + t114) * t139; -t123 * Ifges(6,6) + t125 * t115 + Ifges(5,5); t125 * Ifges(6,6) + t123 * t115 + Ifges(5,6); (t130 * m(7)) + 0.2e1 * pkin(4) * t133 + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t133; t125 * mrSges(6,2) + t123 * t118 + mrSges(5,2); mrSges(6,3) - t135; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
