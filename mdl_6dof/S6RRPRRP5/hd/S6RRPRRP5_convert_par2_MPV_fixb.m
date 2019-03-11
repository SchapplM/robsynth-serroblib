% Return the minimum parameter vector for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRP5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t160 = (m(5) + m(6));
t163 = (pkin(9) ^ 2);
t164 = (pkin(4) ^ 2);
t174 = (t164 * m(6) + Ifges(5,2));
t170 = 2 * pkin(9) * mrSges(5,3) + t174;
t147 = t163 * t160 + Ifges(4,1) + t170;
t165 = (pkin(3) ^ 2);
t151 = t165 * t160 + Ifges(4,2);
t177 = t147 - t151;
t176 = (m(3) * pkin(8));
t175 = (-Ifges(6,2) - Ifges(7,2));
t157 = sin(pkin(11));
t159 = cos(pkin(11));
t173 = t157 * t159;
t152 = t157 ^ 2;
t154 = t159 ^ 2;
t172 = t154 - t152;
t171 = 2 * pkin(10) * mrSges(6,3) - t175;
t169 = pkin(10) * m(6) + mrSges(6,3);
t167 = pkin(9) * t160 + mrSges(5,3);
t148 = t167 * pkin(3) + Ifges(4,4);
t168 = t148 * t173;
t149 = mrSges(4,2) - t167;
t150 = pkin(3) * t160 + mrSges(4,1);
t166 = -t157 * t149 + t159 * t150;
t162 = pkin(10) ^ 2;
t158 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (0.2e1 * t168 + t152 * t147 + t154 * t151 + Ifges(3,2) + ((2 * mrSges(3,3) + t176) * pkin(8))) * t158 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t176) * t158; t172 * t177 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t168; t148 * t172 + t173 * t177 + Ifges(3,4); Ifges(4,5) * t159 - Ifges(4,6) * t157 + Ifges(3,5); Ifges(4,5) * t157 + Ifges(4,6) * t159 + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + ((t163 + t165) * t160) + 0.2e1 * t166 * pkin(2) + t170; mrSges(3,1) + t166; t149 * t159 + t150 * t157 + mrSges(3,2); mrSges(4,3); m(4) + t160; m(6) * t162 + Ifges(5,1) + t171 - t174; pkin(4) * t169 + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t162 + t164) * m(6) + t171; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t169; Ifges(6,1) + Ifges(7,1) + t175; Ifges(6,4) + Ifges(7,4); Ifges(6,5) + Ifges(7,5); Ifges(6,6) + Ifges(7,6); 2 * mrSges(7,1) * pkin(5) + Ifges(6,3) + Ifges(7,3); mrSges(6,1) + mrSges(7,1); mrSges(6,2) + mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
