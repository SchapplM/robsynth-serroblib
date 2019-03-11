% Return the minimum parameter vector for
% S6RRPRRR5
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
% MPV [33x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRR5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR5_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR5_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t170 = (m(6) + m(7));
t164 = (m(5) + t170);
t173 = (pkin(9) ^ 2);
t175 = (pkin(4) ^ 2);
t187 = (t175 * t170 + Ifges(5,2));
t183 = 2 * pkin(9) * mrSges(5,3) + t187;
t151 = t164 * t173 + Ifges(4,1) + t183;
t176 = (pkin(3) ^ 2);
t156 = t164 * t176 + Ifges(4,2);
t189 = t151 - t156;
t188 = (m(3) * pkin(8));
t171 = (pkin(11) ^ 2);
t174 = (pkin(5) ^ 2);
t186 = (Ifges(6,2) + (t171 + t174) * m(7));
t167 = sin(pkin(12));
t169 = cos(pkin(12));
t185 = t167 * t169;
t161 = t167 ^ 2;
t163 = t169 ^ 2;
t184 = t163 - t161;
t182 = -pkin(11) * m(7) - mrSges(7,3);
t180 = pkin(9) * t164 + mrSges(5,3);
t152 = t180 * pkin(3) + Ifges(4,4);
t181 = t152 * t185;
t158 = (mrSges(6,3) - t182);
t179 = pkin(10) * t170 + t158;
t178 = 2 * pkin(11) * mrSges(7,3) + 2 * pkin(10) * t158 + Ifges(7,2) + t186;
t153 = mrSges(4,2) - t180;
t155 = pkin(3) * t164 + mrSges(4,1);
t177 = -t153 * t167 + t155 * t169;
t172 = pkin(10) ^ 2;
t168 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (0.2e1 * t181 + t161 * t151 + t163 * t156 + Ifges(3,2) + ((2 * mrSges(3,3) + t188) * pkin(8))) * t168 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t188) * t168; t189 * t184 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t181; t184 * t152 + t189 * t185 + Ifges(3,4); Ifges(4,5) * t169 - Ifges(4,6) * t167 + Ifges(3,5); Ifges(4,5) * t167 + Ifges(4,6) * t169 + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + ((t173 + t176) * t164) + 0.2e1 * t177 * pkin(2) + t183; mrSges(3,1) + t177; t153 * t169 + t155 * t167 + mrSges(3,2); mrSges(4,3); m(4) + t164; t170 * t172 + Ifges(5,1) + t178 - t187; t179 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t172 + t175) * t170 + t178; pkin(4) * t170 + mrSges(5,1); mrSges(5,2) - t179; m(7) * t171 + Ifges(6,1) - t186; Ifges(6,4); t182 * pkin(5) + Ifges(6,5); Ifges(6,6); m(7) * t174 + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
