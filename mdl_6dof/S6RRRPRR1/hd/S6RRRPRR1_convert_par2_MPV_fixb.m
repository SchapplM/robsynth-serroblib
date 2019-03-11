% Return the minimum parameter vector for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% MPV [33x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t127 = (m(6) + m(7));
t131 = (pkin(9) ^ 2);
t134 = (pkin(4) ^ 2);
t133 = (pkin(5) ^ 2);
t146 = (t133 * m(7) + Ifges(6,2));
t141 = 2 * pkin(9) * mrSges(6,3) + t146;
t116 = Ifges(5,2) + (t131 + t134) * t127 + t141;
t117 = t131 * t127 + Ifges(5,1) + t141;
t147 = -t116 + t117;
t145 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t125 = sin(pkin(11));
t126 = cos(pkin(11));
t144 = t125 * t126;
t121 = t125 ^ 2;
t122 = t126 ^ 2;
t143 = t122 - t121;
t142 = Ifges(5,4) * t144;
t140 = -pkin(8) * m(4) - mrSges(4,3);
t139 = pkin(10) * m(7) + mrSges(7,3);
t138 = -pkin(9) * t127 - mrSges(6,3);
t137 = (mrSges(3,3) - t140);
t119 = pkin(4) * t127 + mrSges(5,1);
t136 = -t125 * mrSges(5,2) + t126 * t119;
t135 = pkin(2) ^ 2;
t132 = pkin(8) ^ 2;
t130 = pkin(10) ^ 2;
t128 = (m(3) + m(4));
t120 = t132 + t135;
t118 = t138 * pkin(4) + Ifges(5,5);
t1 = [Ifges(2,3) + Ifges(3,2) + Ifges(4,2) + t121 * t117 + 0.2e1 * t142 + t122 * t116 + (2 * pkin(8) * mrSges(4,3)) + t120 * m(4) + (2 * pkin(7) * t137) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * t128); pkin(1) * t128 + mrSges(2,1); -pkin(7) * t128 + mrSges(2,2) - t137; Ifges(3,1) - Ifges(3,2) + (-t120 + t132) * m(4); Ifges(3,4); t140 * pkin(2) + Ifges(3,5); Ifges(3,6); t135 * m(4) + Ifges(3,3); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2); t147 * t143 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t142; t143 * Ifges(5,4) + t147 * t144 + Ifges(4,4); -t125 * Ifges(5,6) + t126 * t118 + Ifges(4,5); t126 * Ifges(5,6) + t125 * t118 + Ifges(4,6); 0.2e1 * pkin(3) * t136 + (t134 * t127) + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t136; t126 * mrSges(5,2) + t125 * t119 + mrSges(4,2); mrSges(5,3) - t138; m(5) + t127; m(7) * t130 + Ifges(6,1) + t145 - t146; t139 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t130 + t133) * m(7) + t145; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t139; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
