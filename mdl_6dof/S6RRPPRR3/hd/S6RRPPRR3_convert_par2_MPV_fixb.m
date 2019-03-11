% Return the minimum parameter vector for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% MPV [30x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPPRR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t174 = (m(6) + m(7));
t177 = (pkin(9) ^ 2);
t179 = (pkin(4) ^ 2);
t178 = (pkin(5) ^ 2);
t190 = (t178 * m(7) + Ifges(6,2));
t184 = 2 * pkin(9) * mrSges(6,3) + t190;
t158 = Ifges(5,2) + (t177 + t179) * t174 + t184;
t159 = t177 * t174 + Ifges(5,1) + t184;
t169 = sin(pkin(12));
t162 = t169 ^ 2;
t172 = cos(pkin(12));
t165 = t172 ^ 2;
t188 = t169 * t172;
t185 = Ifges(5,4) * t188;
t155 = t162 * t158 + t165 * t159 + Ifges(4,1) - 0.2e1 * t185;
t161 = t179 * t174 + Ifges(4,2) + Ifges(5,3);
t192 = t155 - t161;
t191 = (m(3) * pkin(8));
t189 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t170 = sin(pkin(11));
t173 = cos(pkin(11));
t187 = t170 * t173;
t163 = t170 ^ 2;
t166 = t173 ^ 2;
t186 = t166 - t163;
t183 = pkin(10) * m(7) + mrSges(7,3);
t181 = -pkin(9) * t174 - mrSges(6,3);
t160 = t181 * pkin(4) + Ifges(5,5);
t157 = t169 * Ifges(5,6) - t172 * t160 + Ifges(4,4);
t182 = t157 * t187;
t180 = t173 * mrSges(4,1) - t170 * mrSges(4,2);
t176 = pkin(10) ^ 2;
t171 = sin(pkin(6));
t156 = -t172 * Ifges(5,6) - t169 * t160 + Ifges(4,6);
t154 = Ifges(4,5) + (t165 - t162) * Ifges(5,4) + (-t158 + t159) * t188;
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (0.2e1 * t182 + t163 * t155 + t166 * t161 + Ifges(3,2) + ((2 * mrSges(3,3) + t191) * pkin(8))) * t171 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t191) * t171; t192 * t186 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t182; t186 * t157 + t192 * t187 + Ifges(3,4); t173 * t154 - t170 * t156 + Ifges(3,5); t170 * t154 + t173 * t156 + Ifges(3,6); 0.2e1 * pkin(2) * t180 + t165 * t158 + t162 * t159 + Ifges(3,3) + Ifges(4,3) + 0.2e1 * t185; mrSges(3,1) + t180; t170 * mrSges(4,1) + t173 * mrSges(4,2) + mrSges(3,2); mrSges(4,3); m(4); pkin(4) * t174 + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t181; m(5) + t174; m(7) * t176 + Ifges(6,1) + t189 - t190; t183 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t176 + t178) * m(7) + t189; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t183; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
