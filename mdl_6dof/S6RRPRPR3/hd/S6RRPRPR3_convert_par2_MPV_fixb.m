% Return the minimum parameter vector for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% MPV [28x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRPR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t144 = pkin(8) ^ 2;
t143 = (pkin(9) ^ 2);
t145 = (pkin(5) ^ 2);
t158 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t124 = Ifges(6,2) + (t143 + t145) * m(7) + t158;
t126 = m(7) * t143 + Ifges(6,1) + t158;
t139 = sin(pkin(11));
t133 = t139 ^ 2;
t141 = cos(pkin(11));
t135 = t141 ^ 2;
t153 = t135 * t124 + t133 * t126 + Ifges(5,2);
t157 = t139 * t141;
t154 = Ifges(6,4) * t157;
t149 = (2 * pkin(8) * mrSges(5,3)) + t153 + 0.2e1 * t154;
t121 = t144 * m(5) + Ifges(4,1) + t149;
t146 = pkin(3) ^ 2;
t131 = t146 * m(5) + Ifges(4,2);
t159 = t121 - t131;
t140 = sin(pkin(10));
t142 = cos(pkin(10));
t156 = t140 * t142;
t134 = t140 ^ 2;
t136 = t142 ^ 2;
t155 = t136 - t134;
t152 = pkin(8) * m(5) + mrSges(5,3);
t151 = -pkin(9) * m(7) - mrSges(7,3);
t128 = t152 * pkin(3) + Ifges(4,4);
t150 = t128 * t156;
t130 = m(7) * pkin(5) + mrSges(6,1);
t148 = -t139 * mrSges(6,2) + t141 * t130;
t129 = mrSges(4,2) - t152;
t132 = m(5) * pkin(3) + mrSges(4,1);
t147 = -t140 * t129 + t142 * t132;
t127 = t151 * pkin(5) + Ifges(6,5);
t1 = [Ifges(2,3) + Ifges(3,2) + t134 * t121 + 0.2e1 * t150 + t136 * t131 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t159 * t155 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t150; t155 * t128 + t159 * t156 + Ifges(3,4); t142 * Ifges(4,5) - t140 * Ifges(4,6) + Ifges(3,5); t140 * Ifges(4,5) + t142 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + (t144 + t146) * m(5) + 0.2e1 * t147 * pkin(2) + t149; mrSges(3,1) + t147; t142 * t129 + t140 * t132 + mrSges(3,2); mrSges(4,3); m(4) + m(5); t133 * t124 + t135 * t126 + Ifges(5,1) - t153 - 0.4e1 * t154; Ifges(5,4) + (t135 - t133) * Ifges(6,4) + (-t124 + t126) * t157; -t139 * Ifges(6,6) + t141 * t127 + Ifges(5,5); t141 * Ifges(6,6) + t139 * t127 + Ifges(5,6); (t145 * m(7)) + 0.2e1 * pkin(4) * t148 + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t148; t141 * mrSges(6,2) + t139 * t130 + mrSges(5,2); mrSges(6,3) - t151; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
