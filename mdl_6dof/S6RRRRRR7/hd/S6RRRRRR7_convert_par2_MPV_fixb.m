% Return the minimum parameter vector for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% MPV [38x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t180 = -pkin(12) * m(7) - mrSges(7,3);
t149 = (mrSges(6,3) - t180);
t160 = (m(6) + m(7));
t179 = -pkin(11) * t160 - t149;
t156 = (m(5) + t160);
t152 = (m(4) + t156);
t150 = m(3) + t152;
t177 = (t150 * pkin(8));
t168 = (pkin(3) ^ 2);
t176 = (t168 * t156 + Ifges(4,2));
t163 = (pkin(11) ^ 2);
t167 = (pkin(4) ^ 2);
t175 = (Ifges(5,2) + (t163 + t167) * t160);
t162 = (pkin(12) ^ 2);
t166 = (pkin(5) ^ 2);
t174 = (Ifges(6,2) + (t162 + t166) * m(7));
t173 = 2 * pkin(9) * mrSges(4,3) + t176;
t172 = pkin(9) * t152 + mrSges(4,3);
t145 = (mrSges(5,3) - t179);
t171 = pkin(10) * t156 + t145;
t170 = 2 * pkin(12) * mrSges(7,3) + 2 * pkin(10) * t145 + 2 * pkin(11) * t149 + Ifges(7,2) + t174 + t175;
t169 = (pkin(2) ^ 2);
t165 = pkin(9) ^ 2;
t164 = pkin(10) ^ 2;
t159 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * t150 + Ifges(2,3) + (t169 * t152 + Ifges(3,2) + (2 * mrSges(3,3) + t177) * pkin(8)) * t159 ^ 2; pkin(1) * t150 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t177) * t159; Ifges(3,1) - Ifges(3,2) + (t165 - t169) * t152 + t173; t172 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t165 + t169) * t152 + t173; pkin(2) * t152 + mrSges(3,1); mrSges(3,2) - t172; t164 * t156 + Ifges(4,1) + t170 - t176; t171 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t164 + t168) * t156 + t170; pkin(3) * t156 + mrSges(4,1); mrSges(4,2) - t171; t163 * t160 + Ifges(5,1) - t175; Ifges(5,4); t179 * pkin(4) + Ifges(5,5); Ifges(5,6); t167 * t160 + Ifges(5,3); pkin(4) * t160 + mrSges(5,1); mrSges(5,2); m(7) * t162 + Ifges(6,1) - t174; Ifges(6,4); t180 * pkin(5) + Ifges(6,5); Ifges(6,6); t166 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
