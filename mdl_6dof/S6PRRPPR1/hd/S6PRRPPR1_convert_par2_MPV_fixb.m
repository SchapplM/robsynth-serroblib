% Return the minimum parameter vector for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% MPV [24x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRPPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t153 = (pkin(9) ^ 2);
t154 = (pkin(5) ^ 2);
t162 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t140 = Ifges(6,2) + (t153 + t154) * m(7) + t162;
t141 = m(7) * t153 + Ifges(6,1) + t162;
t149 = sin(pkin(12));
t144 = t149 ^ 2;
t151 = cos(pkin(12));
t146 = t151 ^ 2;
t161 = t149 * t151;
t158 = Ifges(6,4) * t161;
t137 = t140 * t144 + t141 * t146 + Ifges(5,1) - 0.2e1 * t158;
t143 = m(7) * t154 + Ifges(5,2) + Ifges(6,3);
t163 = t137 - t143;
t150 = sin(pkin(11));
t152 = cos(pkin(11));
t160 = t150 * t152;
t145 = t150 ^ 2;
t147 = t152 ^ 2;
t159 = t147 - t145;
t157 = -m(7) * pkin(9) - mrSges(7,3);
t142 = pkin(5) * t157 + Ifges(6,5);
t139 = Ifges(6,6) * t149 - t142 * t151 + Ifges(5,4);
t156 = t139 * t160;
t155 = mrSges(5,1) * t152 - mrSges(5,2) * t150;
t138 = -Ifges(6,6) * t151 - t142 * t149 + Ifges(5,6);
t136 = Ifges(5,5) + (t146 - t144) * Ifges(6,4) + (-t140 + t141) * t161;
t1 = [m(2) + m(3) + m(4); Ifges(3,3) + Ifges(4,2) + t145 * t137 + 0.2e1 * t156 + t147 * t143 + (2 * pkin(8) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(8) ^ 2) * m(4)); m(4) * pkin(2) + mrSges(3,1); -m(4) * pkin(8) + mrSges(3,2) - mrSges(4,3); t159 * t163 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t156; t139 * t159 + t160 * t163 + Ifges(4,4); t136 * t152 - t138 * t150 + Ifges(4,5); t136 * t150 + t138 * t152 + Ifges(4,6); 0.2e1 * pkin(3) * t155 + t146 * t140 + t144 * t141 + Ifges(4,3) + Ifges(5,3) + 0.2e1 * t158; mrSges(4,1) + t155; mrSges(5,1) * t150 + mrSges(5,2) * t152 + mrSges(4,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t157; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
