% Return the minimum parameter vector for
% S6RRRRPR9
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
% MPV [35x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRPR9_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t165 = (m(4) + m(5));
t159 = m(3) + t165;
t182 = (pkin(8) * t159);
t168 = (pkin(10) ^ 2);
t171 = (pkin(3) ^ 2);
t181 = (Ifges(4,2) + (t168 + t171) * m(5));
t180 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t162 = sin(pkin(12));
t164 = cos(pkin(12));
t179 = t162 * t164;
t170 = (pkin(5) ^ 2);
t178 = (t170 * m(7) + Ifges(5,2) + Ifges(6,3));
t177 = Ifges(6,4) * t179;
t176 = -pkin(10) * m(5) - mrSges(5,3);
t175 = -pkin(11) * m(7) - mrSges(7,3);
t154 = (mrSges(4,3) - t176);
t174 = pkin(9) * t165 + t154;
t173 = 2 * pkin(10) * mrSges(5,3) + 2 * pkin(9) * t154 + t178 + t181;
t172 = (pkin(2) ^ 2);
t169 = pkin(9) ^ 2;
t167 = pkin(11) ^ 2;
t163 = sin(pkin(6));
t158 = t164 ^ 2;
t156 = t162 ^ 2;
t151 = t175 * pkin(5) + Ifges(6,5);
t150 = m(7) * t167 + Ifges(6,1) + t180;
t149 = Ifges(6,2) + (t167 + t170) * m(7) + t180;
t1 = [pkin(1) ^ 2 * t159 + Ifges(2,3) + (t172 * t165 + Ifges(3,2) + (2 * mrSges(3,3) + t182) * pkin(8)) * t163 ^ 2; pkin(1) * t159 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t182) * t163; Ifges(3,1) - Ifges(3,2) + (t169 - t172) * t165 + t173; t174 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t169 + t172) * t165 + t173; pkin(2) * t165 + mrSges(3,1); mrSges(3,2) - t174; m(5) * t168 + Ifges(4,1) - t181; Ifges(4,4); t176 * pkin(3) + Ifges(4,5); Ifges(4,6); m(5) * t171 + Ifges(4,3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); t149 * t156 + t150 * t158 + Ifges(5,1) - 0.2e1 * t177 - t178; Ifges(6,6) * t162 - t151 * t164 + Ifges(5,4); Ifges(5,5) + (t158 - t156) * Ifges(6,4) + (-t149 + t150) * t179; -Ifges(6,6) * t164 - t151 * t162 + Ifges(5,6); t149 * t158 + t150 * t156 + Ifges(5,3) + 0.2e1 * t177; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t175; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
