% Return the minimum parameter vector for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRR8_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR8_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t168 = (m(6) + m(7));
t169 = m(3) + m(4);
t186 = (t169 * pkin(8));
t170 = (pkin(11) ^ 2);
t173 = (pkin(5) ^ 2);
t185 = (Ifges(6,2) + (t170 + t173) * m(7));
t165 = sin(pkin(12));
t167 = cos(pkin(12));
t184 = t165 * t167;
t171 = (pkin(10) ^ 2);
t181 = -pkin(11) * m(7) - mrSges(7,3);
t158 = (mrSges(6,3) - t181);
t178 = 2 * pkin(11) * mrSges(7,3) + 2 * pkin(10) * t158 + Ifges(7,2) + t185;
t149 = t171 * t168 + Ifges(5,1) + t178;
t174 = (pkin(4) ^ 2);
t156 = t174 * t168 + Ifges(5,2);
t160 = t165 ^ 2;
t162 = t167 ^ 2;
t183 = t160 * t149 + t162 * t156 + Ifges(4,2);
t182 = pkin(9) * m(4) + mrSges(4,3);
t179 = pkin(10) * t168 + t158;
t151 = t179 * pkin(4) + Ifges(5,4);
t180 = t151 * t184;
t177 = (2 * pkin(9) * mrSges(4,3)) + 0.2e1 * t180 + t183;
t152 = mrSges(5,2) - t179;
t155 = pkin(4) * t168 + mrSges(5,1);
t176 = -t165 * t152 + t167 * t155;
t175 = (pkin(2) ^ 2);
t172 = pkin(9) ^ 2;
t166 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * t169 + Ifges(2,3) + (t175 * m(4) + Ifges(3,2) + (2 * mrSges(3,3) + t186) * pkin(8)) * t166 ^ 2; pkin(1) * t169 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t186) * t166; Ifges(3,1) - Ifges(3,2) + ((t172 - t175) * m(4)) + t177; t182 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t172 + t175) * m(4)) + t177; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t182; t162 * t149 + t160 * t156 + Ifges(4,1) - 0.4e1 * t180 - t183; Ifges(4,4) + (t162 - t160) * t151 + (t149 - t156) * t184; t167 * Ifges(5,5) - t165 * Ifges(5,6) + Ifges(4,5); t165 * Ifges(5,5) + t167 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t171 + t174) * t168) + 0.2e1 * t176 * pkin(3) + t178; mrSges(4,1) + t176; t167 * t152 + t165 * t155 + mrSges(4,2); mrSges(5,3); m(5) + t168; m(7) * t170 + Ifges(6,1) - t185; Ifges(6,4); pkin(5) * t181 + Ifges(6,5); Ifges(6,6); t173 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
