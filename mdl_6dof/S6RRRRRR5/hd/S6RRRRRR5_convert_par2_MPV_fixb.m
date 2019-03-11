% Return the minimum parameter vector for
% S6RRRRRR5
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRR5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t165 = (m(6) + m(7));
t185 = -pkin(11) * t165 - mrSges(6,3);
t154 = (mrSges(5,3) - t185);
t161 = (m(5) + t165);
t184 = -pkin(10) * t161 - t154;
t156 = (m(4) + t161);
t155 = m(3) + t156;
t182 = (t155 * pkin(8));
t169 = (pkin(10) ^ 2);
t173 = (pkin(3) ^ 2);
t181 = (Ifges(4,2) + (t169 + t173) * t161);
t168 = (pkin(11) ^ 2);
t172 = (pkin(4) ^ 2);
t180 = (Ifges(5,2) + (t168 + t172) * t165);
t171 = (pkin(5) ^ 2);
t179 = (t171 * m(7) + Ifges(6,2));
t178 = 2 * pkin(12) * mrSges(7,3) + Ifges(7,2);
t177 = pkin(12) * m(7) + mrSges(7,3);
t150 = (mrSges(4,3) - t184);
t176 = pkin(9) * t156 + t150;
t175 = 2 * pkin(11) * mrSges(6,3) + 2 * pkin(9) * t150 + 2 * pkin(10) * t154 + t179 + t180 + t181;
t174 = (pkin(2) ^ 2);
t170 = pkin(9) ^ 2;
t167 = pkin(12) ^ 2;
t164 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * t155 + Ifges(2,3) + (t174 * t156 + Ifges(3,2) + (2 * mrSges(3,3) + t182) * pkin(8)) * t164 ^ 2; pkin(1) * t155 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t182) * t164; Ifges(3,1) - Ifges(3,2) + (t170 - t174) * t156 + t175; t176 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t170 + t174) * t156 + t175; pkin(2) * t156 + mrSges(3,1); mrSges(3,2) - t176; t169 * t161 + Ifges(4,1) - t181; Ifges(4,4); t184 * pkin(3) + Ifges(4,5); Ifges(4,6); t173 * t161 + Ifges(4,3); pkin(3) * t161 + mrSges(4,1); mrSges(4,2); t168 * t165 + Ifges(5,1) - t180; Ifges(5,4); t185 * pkin(4) + Ifges(5,5); Ifges(5,6); t172 * t165 + Ifges(5,3); pkin(4) * t165 + mrSges(5,1); mrSges(5,2); m(7) * t167 + Ifges(6,1) + t178 - t179; t177 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t167 + t171) * m(7) + t178; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t177; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
