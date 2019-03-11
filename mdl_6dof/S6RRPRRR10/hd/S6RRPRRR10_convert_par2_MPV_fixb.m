% Return the minimum parameter vector for
% S6RRPRRR10
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t181 = (m(3) * pkin(8));
t164 = (m(6) + m(7));
t169 = (pkin(4) ^ 2);
t180 = (t169 * t164 + Ifges(5,2));
t165 = (pkin(11) ^ 2);
t168 = (pkin(5) ^ 2);
t179 = (Ifges(6,2) + (t165 + t168) * m(7));
t161 = sin(pkin(12));
t163 = cos(pkin(12));
t178 = t161 * t163;
t158 = (m(5) + t164);
t177 = Ifges(4,4) * t178;
t176 = 2 * pkin(9) * mrSges(5,3) + t180;
t175 = -pkin(11) * m(7) - mrSges(7,3);
t174 = -pkin(9) * t158 - mrSges(5,3);
t152 = (mrSges(6,3) - t175);
t173 = pkin(10) * t164 + t152;
t170 = (pkin(3) ^ 2);
t172 = t170 * t158 + Ifges(3,2) + Ifges(4,3);
t171 = 2 * pkin(11) * mrSges(7,3) + 2 * pkin(10) * t152 + Ifges(7,2) + t179;
t167 = pkin(9) ^ 2;
t166 = pkin(10) ^ 2;
t162 = sin(pkin(6));
t157 = t163 ^ 2;
t155 = t161 ^ 2;
t149 = t174 * pkin(3) + Ifges(4,5);
t148 = t167 * t158 + Ifges(4,1) + t176;
t147 = Ifges(4,2) + (t167 + t170) * t158 + t176;
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + ((2 * mrSges(3,3) + t181) * pkin(8) + t172) * t162 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t181) * t162; t155 * t147 + t157 * t148 + Ifges(3,1) - t172 - 0.2e1 * t177; t161 * Ifges(4,6) - t163 * t149 + Ifges(3,4); Ifges(3,5) + (t157 - t155) * Ifges(4,4) + (-t147 + t148) * t178; -t163 * Ifges(4,6) - t161 * t149 + Ifges(3,6); t157 * t147 + t155 * t148 + Ifges(3,3) + 0.2e1 * t177; mrSges(3,1); mrSges(3,2); pkin(3) * t158 + mrSges(4,1); mrSges(4,2); mrSges(4,3) - t174; m(4) + t158; t166 * t164 + Ifges(5,1) + t171 - t180; t173 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t166 + t169) * t164 + t171; pkin(4) * t164 + mrSges(5,1); mrSges(5,2) - t173; m(7) * t165 + Ifges(6,1) - t179; Ifges(6,4); t175 * pkin(5) + Ifges(6,5); Ifges(6,6); t168 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
