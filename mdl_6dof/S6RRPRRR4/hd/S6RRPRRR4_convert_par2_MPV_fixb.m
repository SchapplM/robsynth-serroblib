% Return the minimum parameter vector for
% S6RRPRRR4
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t172 = (m(6) + m(7));
t166 = (m(5) + t172);
t176 = (pkin(9) ^ 2);
t183 = -pkin(10) * t172 - mrSges(6,3);
t159 = (mrSges(5,3) - t183);
t177 = (pkin(5) ^ 2);
t189 = (t177 * m(7) + Ifges(6,2));
t175 = (pkin(10) ^ 2);
t178 = (pkin(4) ^ 2);
t190 = (Ifges(5,2) + (t175 + t178) * t172);
t180 = 2 * pkin(10) * mrSges(6,3) + 2 * pkin(9) * t159 + t189 + t190;
t154 = t176 * t166 + Ifges(4,1) + t180;
t179 = (pkin(3) ^ 2);
t161 = t179 * t166 + Ifges(4,2);
t192 = t154 - t161;
t191 = (m(3) * pkin(8));
t188 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t169 = sin(pkin(12));
t171 = cos(pkin(12));
t187 = t169 * t171;
t163 = t169 ^ 2;
t165 = t171 ^ 2;
t186 = t165 - t163;
t185 = pkin(11) * m(7) + mrSges(7,3);
t182 = pkin(9) * t166 + t159;
t155 = pkin(3) * t182 + Ifges(4,4);
t184 = t155 * t187;
t156 = mrSges(4,2) - t182;
t160 = pkin(3) * t166 + mrSges(4,1);
t181 = -t169 * t156 + t171 * t160;
t174 = pkin(11) ^ 2;
t170 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (0.2e1 * t184 + t163 * t154 + t165 * t161 + Ifges(3,2) + ((2 * mrSges(3,3) + t191) * pkin(8))) * t170 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t191) * t170; t192 * t186 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t184; t186 * t155 + t192 * t187 + Ifges(3,4); t171 * Ifges(4,5) - t169 * Ifges(4,6) + Ifges(3,5); t169 * Ifges(4,5) + t171 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + ((t176 + t179) * t166) + 0.2e1 * t181 * pkin(2) + t180; mrSges(3,1) + t181; t171 * t156 + t169 * t160 + mrSges(3,2); mrSges(4,3); m(4) + t166; t175 * t172 + Ifges(5,1) - t190; Ifges(5,4); pkin(4) * t183 + Ifges(5,5); Ifges(5,6); t178 * t172 + Ifges(5,3); pkin(4) * t172 + mrSges(5,1); mrSges(5,2); m(7) * t174 + Ifges(6,1) + t188 - t189; pkin(5) * t185 + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t174 + t177) * m(7) + t188; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t185; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
