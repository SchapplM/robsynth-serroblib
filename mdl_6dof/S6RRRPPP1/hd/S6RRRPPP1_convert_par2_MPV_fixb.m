% Return the minimum parameter vector for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% MPV [29x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPPP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t166 = Ifges(5,2) + Ifges(6,3) + Ifges(7,2);
t169 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t175 = sin(pkin(10));
t170 = t175 ^ 2;
t177 = cos(pkin(10));
t172 = t177 ^ 2;
t199 = t166 * t172 + t169 * t170;
t168 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t198 = (t172 - t170) * t168;
t165 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t176 = sin(pkin(6));
t197 = t165 * t176;
t167 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t195 = t167 * t175;
t178 = cos(pkin(6));
t193 = t176 * t178;
t192 = t177 * t178;
t184 = t192 * t197;
t185 = t193 * t195;
t191 = 0.2e1 * t185 + 0.2e1 * t184;
t190 = pkin(9) * m(4) + mrSges(4,3);
t189 = t177 * t175 * t168;
t188 = (-t166 + t169) * t177;
t187 = 0.2e1 * t189;
t164 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t171 = t176 ^ 2;
t173 = t178 ^ 2;
t186 = t171 * t164 + t199 * t173 + Ifges(4,2);
t183 = (2 * pkin(9) * mrSges(4,3)) + t173 * t187 - 0.2e1 * t184 - 0.2e1 * t185 + t186;
t182 = t187 + t199;
t181 = (pkin(2) ^ 2);
t180 = pkin(9) ^ 2;
t179 = (m(3) + m(4));
t1 = [Ifges(2,3) + t181 * m(4) + Ifges(3,2) + 2 * pkin(8) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(8) ^ 2) * t179; pkin(1) * t179 + mrSges(2,1); -pkin(8) * t179 + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + ((t180 - t181) * m(4)) + t183; t190 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t180 + t181) * m(4)) + t183; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t190; t170 * t166 + t172 * t169 + Ifges(4,1) + (-0.2e1 * t173 - 0.2e1) * t189 - t186 + t191; t178 * t198 - t177 * t176 * t167 + Ifges(4,4) + (t178 * t188 + t197) * t175; t176 * t198 + t167 * t192 + Ifges(4,5) + (-t165 * t178 + t176 * t188) * t175; Ifges(4,6) + (t165 * t177 + t195) * (t173 - t171) + (-t164 + t182) * t193; t173 * t164 + t182 * t171 + Ifges(4,3) + t191; mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
