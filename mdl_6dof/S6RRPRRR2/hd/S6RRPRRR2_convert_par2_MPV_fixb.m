% Return the minimum parameter vector for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t138 = (m(6) + m(7));
t133 = (m(5) + t138);
t141 = (pkin(8) ^ 2);
t144 = (pkin(3) ^ 2);
t143 = (pkin(4) ^ 2);
t155 = (t143 * t138 + Ifges(5,2));
t150 = 2 * pkin(8) * mrSges(5,3) + t155;
t122 = Ifges(4,2) + (t141 + t144) * t133 + t150;
t123 = t141 * t133 + Ifges(4,1) + t150;
t156 = -t122 + t123;
t139 = (pkin(10) ^ 2);
t142 = (pkin(5) ^ 2);
t154 = (Ifges(6,2) + (t139 + t142) * m(7));
t136 = sin(pkin(11));
t137 = cos(pkin(11));
t153 = t136 * t137;
t131 = t136 ^ 2;
t132 = t137 ^ 2;
t152 = t132 - t131;
t151 = Ifges(4,4) * t153;
t149 = -pkin(10) * m(7) - mrSges(7,3);
t148 = -pkin(8) * t133 - mrSges(5,3);
t128 = (mrSges(6,3) - t149);
t147 = pkin(9) * t138 + t128;
t146 = 2 * pkin(10) * mrSges(7,3) + 2 * pkin(9) * t128 + Ifges(7,2) + t154;
t126 = pkin(3) * t133 + mrSges(4,1);
t145 = -t136 * mrSges(4,2) + t137 * t126;
t140 = pkin(9) ^ 2;
t124 = t148 * pkin(3) + Ifges(4,5);
t1 = [Ifges(2,3) + Ifges(3,2) + t131 * t123 + 0.2e1 * t151 + t132 * t122 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t156 * t152 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t151; t152 * Ifges(4,4) + t156 * t153 + Ifges(3,4); -t136 * Ifges(4,6) + t137 * t124 + Ifges(3,5); t137 * Ifges(4,6) + t136 * t124 + Ifges(3,6); 0.2e1 * pkin(2) * t145 + (t144 * t133) + Ifges(3,3) + Ifges(4,3); mrSges(3,1) + t145; t137 * mrSges(4,2) + t136 * t126 + mrSges(3,2); mrSges(4,3) - t148; m(4) + t133; t140 * t138 + Ifges(5,1) + t146 - t155; t147 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t140 + t143) * t138 + t146; pkin(4) * t138 + mrSges(5,1); mrSges(5,2) - t147; m(7) * t139 + Ifges(6,1) - t154; Ifges(6,4); t149 * pkin(5) + Ifges(6,5); Ifges(6,6); t142 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
