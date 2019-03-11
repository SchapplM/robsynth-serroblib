% Return the minimum parameter vector for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% MPV [30x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPPRR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t156 = -pkin(9) * m(7) - mrSges(7,3);
t139 = (m(6) + m(7));
t142 = (pkin(8) ^ 2);
t144 = (pkin(4) ^ 2);
t128 = (mrSges(6,3) - t156);
t141 = (pkin(9) ^ 2);
t143 = (pkin(5) ^ 2);
t152 = (Ifges(6,2) + (t141 + t143) * m(7));
t146 = 2 * pkin(9) * mrSges(7,3) + 2 * pkin(8) * t128 + Ifges(7,2) + t152;
t120 = Ifges(5,2) + (t142 + t144) * t139 + t146;
t123 = t142 * t139 + Ifges(5,1) + t146;
t135 = sin(pkin(11));
t130 = t135 ^ 2;
t137 = cos(pkin(11));
t132 = t137 ^ 2;
t151 = t135 * t137;
t148 = Ifges(5,4) * t151;
t119 = t130 * t120 + t132 * t123 + Ifges(4,1) - 0.2e1 * t148;
t126 = t144 * t139 + Ifges(4,2) + Ifges(5,3);
t155 = t119 - t126;
t154 = -pkin(8) * t139 - t128;
t136 = sin(pkin(10));
t138 = cos(pkin(10));
t150 = t136 * t138;
t131 = t136 ^ 2;
t133 = t138 ^ 2;
t149 = t133 - t131;
t124 = t154 * pkin(4) + Ifges(5,5);
t122 = t135 * Ifges(5,6) - t137 * t124 + Ifges(4,4);
t147 = t122 * t150;
t145 = t138 * mrSges(4,1) - t136 * mrSges(4,2);
t121 = -t137 * Ifges(5,6) - t135 * t124 + Ifges(4,6);
t118 = Ifges(4,5) + (t132 - t130) * Ifges(5,4) + (-t120 + t123) * t151;
t1 = [Ifges(2,3) + Ifges(3,2) + t131 * t119 + 0.2e1 * t147 + t133 * t126 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t155 * t149 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t147; t149 * t122 + t155 * t150 + Ifges(3,4); t138 * t118 - t136 * t121 + Ifges(3,5); t136 * t118 + t138 * t121 + Ifges(3,6); 0.2e1 * pkin(2) * t145 + t132 * t120 + t130 * t123 + Ifges(3,3) + Ifges(4,3) + 0.2e1 * t148; mrSges(3,1) + t145; t136 * mrSges(4,1) + t138 * mrSges(4,2) + mrSges(3,2); mrSges(4,3); m(4); pkin(4) * t139 + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t154; m(5) + t139; m(7) * t141 + Ifges(6,1) - t152; Ifges(6,4); t156 * pkin(5) + Ifges(6,5); Ifges(6,6); t143 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
