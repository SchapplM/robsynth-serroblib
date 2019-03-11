% Return the minimum parameter vector for
% S6RRRPRR12
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
% MPV [35x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRR12_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t179 = -pkin(11) * m(7) - mrSges(7,3);
t149 = (mrSges(6,3) - t179);
t160 = (m(6) + m(7));
t178 = -pkin(10) * t160 - t149;
t161 = m(3) + m(4);
t176 = (t161 * pkin(8));
t163 = (pkin(11) ^ 2);
t166 = (pkin(5) ^ 2);
t175 = (Ifges(6,2) + (t163 + t166) * m(7));
t157 = sin(pkin(12));
t159 = cos(pkin(12));
t174 = t157 * t159;
t167 = (pkin(4) ^ 2);
t173 = (t167 * t160 + Ifges(4,2) + Ifges(5,3));
t172 = Ifges(5,4) * t174;
t171 = pkin(9) * m(4) + mrSges(4,3);
t170 = 2 * pkin(9) * mrSges(4,3) + t173;
t169 = 2 * pkin(11) * mrSges(7,3) + 2 * pkin(10) * t149 + Ifges(7,2) + t175;
t168 = (pkin(2) ^ 2);
t165 = pkin(9) ^ 2;
t164 = pkin(10) ^ 2;
t158 = sin(pkin(6));
t154 = t159 ^ 2;
t152 = t157 ^ 2;
t146 = t178 * pkin(4) + Ifges(5,5);
t145 = t164 * t160 + Ifges(5,1) + t169;
t144 = Ifges(5,2) + (t164 + t167) * t160 + t169;
t1 = [pkin(1) ^ 2 * t161 + Ifges(2,3) + (t168 * m(4) + Ifges(3,2) + (2 * mrSges(3,3) + t176) * pkin(8)) * t158 ^ 2; pkin(1) * t161 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t176) * t158; Ifges(3,1) - Ifges(3,2) + (t165 - t168) * m(4) + t170; t171 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t165 + t168) * m(4) + t170; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t171; t152 * t144 + t154 * t145 + Ifges(4,1) - 0.2e1 * t172 - t173; t157 * Ifges(5,6) - t159 * t146 + Ifges(4,4); Ifges(4,5) + (t154 - t152) * Ifges(5,4) + (-t144 + t145) * t174; -t159 * Ifges(5,6) - t157 * t146 + Ifges(4,6); t154 * t144 + t152 * t145 + Ifges(4,3) + 0.2e1 * t172; mrSges(4,1); mrSges(4,2); pkin(4) * t160 + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t178; m(5) + t160; m(7) * t163 + Ifges(6,1) - t175; Ifges(6,4); t179 * pkin(5) + Ifges(6,5); Ifges(6,6); t166 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
