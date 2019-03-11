% Return the minimum parameter vector for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPPR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t161 = m(3) + m(4);
t174 = (t161 * pkin(8));
t173 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t156 = sin(pkin(11));
t158 = cos(pkin(11));
t172 = t156 * t158;
t162 = pkin(10) ^ 2;
t164 = (pkin(5) ^ 2);
t146 = Ifges(5,2) + Ifges(6,3) + (t162 + t164) * m(7) + t173;
t150 = t164 * m(7) + Ifges(5,1) + Ifges(6,2);
t151 = t156 ^ 2;
t153 = t158 ^ 2;
t171 = t153 * t146 + t151 * t150 + Ifges(4,2);
t170 = pkin(9) * m(4) + mrSges(4,3);
t169 = -pkin(10) * m(7) - mrSges(7,3);
t160 = Ifges(5,4) + Ifges(6,6);
t168 = t160 * t172;
t167 = (2 * pkin(9) * mrSges(4,3)) + 0.2e1 * t168 + t171;
t166 = t158 * mrSges(5,1) - t156 * mrSges(5,2);
t165 = (pkin(2) ^ 2);
t163 = pkin(9) ^ 2;
t159 = Ifges(5,6) - Ifges(6,5);
t157 = sin(pkin(6));
t149 = t169 * pkin(5) - Ifges(6,4) + Ifges(5,5);
t1 = [pkin(1) ^ 2 * t161 + Ifges(2,3) + (t165 * m(4) + Ifges(3,2) + (2 * mrSges(3,3) + t174) * pkin(8)) * t157 ^ 2; pkin(1) * t161 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t174) * t157; Ifges(3,1) - Ifges(3,2) + ((t163 - t165) * m(4)) + t167; t170 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t163 + t165) * m(4)) + t167; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t170; t151 * t146 + t153 * t150 + Ifges(4,1) - 0.4e1 * t168 - t171; Ifges(4,4) + (t153 - t151) * t160 + (-t146 + t150) * t172; t158 * t149 - t156 * t159 + Ifges(4,5); t156 * t149 + t158 * t159 + Ifges(4,6); (m(7) * t162) + 0.2e1 * pkin(3) * t166 + Ifges(6,1) + Ifges(4,3) + Ifges(5,3) + t173; mrSges(4,1) + t166; t156 * mrSges(5,1) + t158 * mrSges(5,2) + mrSges(4,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) + t169; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
