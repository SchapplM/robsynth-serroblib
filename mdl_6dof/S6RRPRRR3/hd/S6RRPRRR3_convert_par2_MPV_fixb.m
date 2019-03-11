% Return the minimum parameter vector for
% S6RRPRRR3
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
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t162 = -pkin(10) * m(7) - mrSges(7,3);
t143 = (m(6) + m(7));
t139 = (m(5) + t143);
t147 = (pkin(8) ^ 2);
t134 = (mrSges(6,3) - t162);
t160 = -pkin(9) * t143 - t134;
t128 = (mrSges(5,3) - t160);
t145 = (pkin(10) ^ 2);
t148 = (pkin(5) ^ 2);
t157 = (Ifges(6,2) + (t145 + t148) * m(7));
t146 = (pkin(9) ^ 2);
t149 = (pkin(4) ^ 2);
t158 = (Ifges(5,2) + (t146 + t149) * t143);
t151 = 2 * pkin(10) * mrSges(7,3) + 2 * pkin(8) * t128 + 2 * pkin(9) * t134 + Ifges(7,2) + t157 + t158;
t124 = t147 * t139 + Ifges(4,1) + t151;
t150 = (pkin(3) ^ 2);
t132 = t150 * t139 + Ifges(4,2);
t161 = t124 - t132;
t141 = sin(pkin(11));
t142 = cos(pkin(11));
t156 = t141 * t142;
t137 = t141 ^ 2;
t138 = t142 ^ 2;
t155 = t138 - t137;
t153 = pkin(8) * t139 + t128;
t125 = t153 * pkin(3) + Ifges(4,4);
t154 = t125 * t156;
t126 = mrSges(4,2) - t153;
t131 = pkin(3) * t139 + mrSges(4,1);
t152 = -t141 * t126 + t142 * t131;
t1 = [Ifges(2,3) + Ifges(3,2) + t137 * t124 + 0.2e1 * t154 + t138 * t132 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t161 * t155 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t154; t155 * t125 + t161 * t156 + Ifges(3,4); t142 * Ifges(4,5) - t141 * Ifges(4,6) + Ifges(3,5); t141 * Ifges(4,5) + t142 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + ((t147 + t150) * t139) + 0.2e1 * t152 * pkin(2) + t151; mrSges(3,1) + t152; t142 * t126 + t141 * t131 + mrSges(3,2); mrSges(4,3); m(4) + t139; t146 * t143 + Ifges(5,1) - t158; Ifges(5,4); t160 * pkin(4) + Ifges(5,5); Ifges(5,6); t149 * t143 + Ifges(5,3); pkin(4) * t143 + mrSges(5,1); mrSges(5,2); m(7) * t145 + Ifges(6,1) - t157; Ifges(6,4); t162 * pkin(5) + Ifges(6,5); Ifges(6,6); t148 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
