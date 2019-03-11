% Return the minimum parameter vector for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t169 = (m(4) + m(5));
t163 = m(3) + t169;
t186 = (t163 * pkin(8));
t171 = pkin(10) ^ 2;
t174 = pkin(3) ^ 2;
t185 = Ifges(4,2) + (t171 + t174) * m(5);
t184 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t166 = sin(pkin(12));
t168 = cos(pkin(12));
t183 = t166 * t168;
t170 = pkin(11) ^ 2;
t150 = m(7) * t170 + Ifges(6,1) + t184;
t173 = (pkin(5) ^ 2);
t157 = t173 * m(7) + Ifges(6,2);
t160 = t166 ^ 2;
t162 = t168 ^ 2;
t182 = t160 * t150 + t162 * t157 + Ifges(5,2);
t181 = -pkin(10) * m(5) - mrSges(5,3);
t180 = pkin(11) * m(7) + mrSges(7,3);
t152 = t180 * pkin(5) + Ifges(6,4);
t179 = t152 * t183;
t156 = (mrSges(4,3) - t181);
t178 = pkin(9) * t169 + t156;
t155 = mrSges(6,2) - t180;
t158 = m(7) * pkin(5) + mrSges(6,1);
t177 = -t166 * t155 + t168 * t158;
t176 = (2 * pkin(10) * mrSges(5,3)) + (2 * pkin(9) * t156) + 0.2e1 * t179 + t182 + t185;
t175 = (pkin(2) ^ 2);
t172 = pkin(9) ^ 2;
t167 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * t163 + Ifges(2,3) + (t175 * t169 + Ifges(3,2) + (2 * mrSges(3,3) + t186) * pkin(8)) * t167 ^ 2; pkin(1) * t163 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t186) * t167; Ifges(3,1) - Ifges(3,2) + ((t172 - t175) * t169) + t176; t178 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t172 + t175) * t169) + t176; pkin(2) * t169 + mrSges(3,1); mrSges(3,2) - t178; t171 * m(5) + Ifges(4,1) - t185; Ifges(4,4); t181 * pkin(3) + Ifges(4,5); Ifges(4,6); t174 * m(5) + Ifges(4,3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); t162 * t150 + t160 * t157 + Ifges(5,1) - 0.4e1 * t179 - t182; Ifges(5,4) + (t162 - t160) * t152 + (t150 - t157) * t183; t168 * Ifges(6,5) - t166 * Ifges(6,6) + Ifges(5,5); t166 * Ifges(6,5) + t168 * Ifges(6,6) + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t170 + t173) * m(7)) + 0.2e1 * t177 * pkin(4) + t184; mrSges(5,1) + t177; t168 * t155 + t166 * t158 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
