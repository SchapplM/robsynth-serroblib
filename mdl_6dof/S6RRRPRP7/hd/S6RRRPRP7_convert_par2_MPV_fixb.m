% Return the minimum parameter vector for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRP7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t168 = m(3) + m(4);
t182 = (t168 * pkin(8));
t181 = (-Ifges(6,2) - Ifges(7,3));
t165 = sin(pkin(11));
t167 = cos(pkin(11));
t180 = t165 * t167;
t179 = 2 * pkin(10) * mrSges(6,3) - t181;
t169 = pkin(10) ^ 2;
t154 = t169 * m(6) + Ifges(5,1) + t179;
t171 = (pkin(4) ^ 2);
t158 = t171 * m(6) + Ifges(5,2);
t160 = t165 ^ 2;
t162 = t167 ^ 2;
t178 = t160 * t154 + t162 * t158 + Ifges(4,2);
t177 = pkin(9) * m(4) + mrSges(4,3);
t176 = pkin(10) * m(6) + mrSges(6,3);
t156 = t176 * pkin(4) + Ifges(5,4);
t175 = t156 * t180;
t174 = (2 * pkin(9) * mrSges(4,3)) + 0.2e1 * t175 + t178;
t157 = mrSges(5,2) - t176;
t159 = m(6) * pkin(4) + mrSges(5,1);
t173 = -t165 * t157 + t167 * t159;
t172 = (pkin(2) ^ 2);
t170 = pkin(9) ^ 2;
t166 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * t168 + Ifges(2,3) + (t172 * m(4) + Ifges(3,2) + (2 * mrSges(3,3) + t182) * pkin(8)) * t166 ^ 2; pkin(1) * t168 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t182) * t166; Ifges(3,1) - Ifges(3,2) + ((t170 - t172) * m(4)) + t174; t177 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t170 + t172) * m(4)) + t174; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t177; t162 * t154 + t160 * t158 + Ifges(4,1) - 0.4e1 * t175 - t178; Ifges(4,4) + (t162 - t160) * t156 + (t154 - t158) * t180; t167 * Ifges(5,5) - t165 * Ifges(5,6) + Ifges(4,5); t165 * Ifges(5,5) + t167 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t169 + t171) * m(6)) + 0.2e1 * t173 * pkin(3) + t179; mrSges(4,1) + t173; t167 * t157 + t165 * t159 + mrSges(4,2); mrSges(5,3); m(5) + m(6); Ifges(6,1) + Ifges(7,1) + t181; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
