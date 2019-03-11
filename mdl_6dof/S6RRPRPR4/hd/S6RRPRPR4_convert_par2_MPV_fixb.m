% Return the minimum parameter vector for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRPR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t180 = pkin(9) ^ 2;
t188 = pkin(10) * m(7) + mrSges(7,3);
t159 = pkin(5) * t188 + Ifges(6,4);
t174 = sin(pkin(12));
t177 = cos(pkin(12));
t193 = t174 * t177;
t187 = t159 * t193;
t179 = (pkin(10) ^ 2);
t194 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t157 = m(7) * t179 + Ifges(6,1) + t194;
t181 = (pkin(5) ^ 2);
t163 = t181 * m(7) + Ifges(6,2);
t167 = t174 ^ 2;
t170 = t177 ^ 2;
t190 = t167 * t157 + t170 * t163 + Ifges(5,2);
t185 = (2 * pkin(9) * mrSges(5,3)) + 0.2e1 * t187 + t190;
t154 = t180 * m(5) + Ifges(4,1) + t185;
t182 = pkin(3) ^ 2;
t165 = t182 * m(5) + Ifges(4,2);
t196 = t154 - t165;
t195 = (m(3) * pkin(8));
t175 = sin(pkin(11));
t178 = cos(pkin(11));
t192 = t175 * t178;
t168 = t175 ^ 2;
t171 = t178 ^ 2;
t191 = t171 - t168;
t189 = pkin(9) * m(5) + mrSges(5,3);
t160 = pkin(3) * t189 + Ifges(4,4);
t186 = t160 * t192;
t161 = mrSges(6,2) - t188;
t164 = m(7) * pkin(5) + mrSges(6,1);
t184 = -t174 * t161 + t177 * t164;
t162 = mrSges(4,2) - t189;
t166 = m(5) * pkin(3) + mrSges(4,1);
t183 = -t175 * t162 + t178 * t166;
t176 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (0.2e1 * t186 + t168 * t154 + t171 * t165 + Ifges(3,2) + ((2 * mrSges(3,3) + t195) * pkin(8))) * t176 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t195) * t176; t196 * t191 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t186; t191 * t160 + t196 * t192 + Ifges(3,4); t178 * Ifges(4,5) - t175 * Ifges(4,6) + Ifges(3,5); t175 * Ifges(4,5) + t178 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + (t180 + t182) * m(5) + 0.2e1 * t183 * pkin(2) + t185; mrSges(3,1) + t183; t178 * t162 + t175 * t166 + mrSges(3,2); mrSges(4,3); m(4) + m(5); t170 * t157 + t167 * t163 + Ifges(5,1) - 0.4e1 * t187 - t190; Ifges(5,4) + (t170 - t167) * t159 + (t157 - t163) * t193; t177 * Ifges(6,5) - t174 * Ifges(6,6) + Ifges(5,5); t174 * Ifges(6,5) + t177 * Ifges(6,6) + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t179 + t181) * m(7)) + 0.2e1 * t184 * pkin(4) + t194; mrSges(5,1) + t184; t177 * t161 + t174 * t164 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
