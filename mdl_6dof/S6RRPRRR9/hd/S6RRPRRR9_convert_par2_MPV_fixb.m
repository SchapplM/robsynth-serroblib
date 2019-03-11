% Return the minimum parameter vector for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% MPV [35x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRR9_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t164 = (m(6) + m(7));
t183 = -pkin(10) * t164 - mrSges(6,3);
t152 = (mrSges(5,3) - t183);
t158 = (m(5) + t164);
t182 = -pkin(9) * t158 - t152;
t181 = (m(3) * pkin(8));
t167 = (pkin(10) ^ 2);
t170 = (pkin(4) ^ 2);
t179 = (Ifges(5,2) + (t167 + t170) * t164);
t169 = (pkin(5) ^ 2);
t178 = (t169 * m(7) + Ifges(6,2));
t177 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t161 = sin(pkin(12));
t163 = cos(pkin(12));
t176 = t161 * t163;
t175 = Ifges(4,4) * t176;
t174 = pkin(11) * m(7) + mrSges(7,3);
t171 = (pkin(3) ^ 2);
t173 = (t171 * t158 + Ifges(3,2) + Ifges(4,3));
t172 = 2 * pkin(10) * mrSges(6,3) + 2 * pkin(9) * t152 + t178 + t179;
t168 = pkin(9) ^ 2;
t166 = pkin(11) ^ 2;
t162 = sin(pkin(6));
t157 = t163 ^ 2;
t155 = t161 ^ 2;
t149 = pkin(3) * t182 + Ifges(4,5);
t148 = t168 * t158 + Ifges(4,1) + t172;
t147 = Ifges(4,2) + (t168 + t171) * t158 + t172;
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + ((2 * mrSges(3,3) + t181) * pkin(8) + t173) * t162 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t181) * t162; t155 * t147 + t157 * t148 + Ifges(3,1) - t173 - 0.2e1 * t175; t161 * Ifges(4,6) - t163 * t149 + Ifges(3,4); Ifges(3,5) + (t157 - t155) * Ifges(4,4) + (-t147 + t148) * t176; -t163 * Ifges(4,6) - t161 * t149 + Ifges(3,6); t157 * t147 + t155 * t148 + Ifges(3,3) + 0.2e1 * t175; mrSges(3,1); mrSges(3,2); pkin(3) * t158 + mrSges(4,1); mrSges(4,2); mrSges(4,3) - t182; m(4) + t158; t167 * t164 + Ifges(5,1) - t179; Ifges(5,4); pkin(4) * t183 + Ifges(5,5); Ifges(5,6); t170 * t164 + Ifges(5,3); pkin(4) * t164 + mrSges(5,1); mrSges(5,2); m(7) * t166 + Ifges(6,1) + t177 - t178; t174 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t166 + t169) * m(7) + t177; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t174; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
