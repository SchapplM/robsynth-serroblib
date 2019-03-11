% Return the minimum parameter vector for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% MPV [32x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPPRR11_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t168 = (m(3) * pkin(8));
t154 = (m(6) + m(7));
t158 = (pkin(5) ^ 2);
t167 = (t158 * m(7) + Ifges(6,2));
t166 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t151 = sin(pkin(11));
t153 = cos(pkin(11));
t165 = t151 * t153;
t164 = Ifges(5,4) * t165;
t163 = 2 * pkin(9) * mrSges(6,3) + t167;
t162 = pkin(10) * m(7) + mrSges(7,3);
t161 = -pkin(9) * t154 - mrSges(6,3);
t157 = pkin(9) ^ 2;
t159 = (pkin(4) ^ 2);
t142 = Ifges(5,2) + (t157 + t159) * t154 + t163;
t143 = t157 * t154 + Ifges(5,1) + t163;
t146 = t151 ^ 2;
t148 = t153 ^ 2;
t160 = t148 * t142 + t146 * t143 + Ifges(3,2) + Ifges(4,3);
t156 = pkin(10) ^ 2;
t152 = sin(pkin(6));
t145 = -0.2e1 * t164;
t144 = t161 * pkin(4) + Ifges(5,5);
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (0.2e1 * t164 + ((2 * mrSges(3,3) + t168) * pkin(8)) + t160) * t152 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t168) * t152; (t159 * t154) + Ifges(3,1) + Ifges(4,2) + Ifges(5,3) + t145 - t160; -t153 * Ifges(5,6) - t151 * t144 + Ifges(3,4) + Ifges(4,6); -t151 * Ifges(5,6) + t153 * t144 - Ifges(4,4) + Ifges(3,5); Ifges(3,6) - Ifges(4,5) - (t148 - t146) * Ifges(5,4) + (t142 - t143) * t165; t146 * t142 + t148 * t143 + Ifges(4,1) + Ifges(3,3) + t145; mrSges(3,1); mrSges(3,2); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); pkin(4) * t154 + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t161; m(5) + t154; m(7) * t156 + Ifges(6,1) + t166 - t167; t162 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t156 + t158) * m(7) + t166; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t162; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
