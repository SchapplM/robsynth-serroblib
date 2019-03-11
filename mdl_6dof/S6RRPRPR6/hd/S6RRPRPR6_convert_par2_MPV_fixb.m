% Return the minimum parameter vector for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRPR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t168 = (pkin(9) ^ 2);
t179 = (pkin(10) * mrSges(7,3));
t162 = 2 * t179;
t167 = (pkin(10) ^ 2);
t169 = (pkin(5) ^ 2);
t173 = (Ifges(5,2) + Ifges(7,2) + Ifges(6,3) + (t167 + t169) * m(7));
t171 = 2 * pkin(9) * mrSges(5,3) + t162 + t173;
t152 = t168 * m(5) + Ifges(4,1) + t171;
t170 = (pkin(3) ^ 2);
t156 = t170 * m(5) + Ifges(4,2);
t181 = t152 - t156;
t180 = (m(3) * pkin(8));
t164 = sin(pkin(11));
t166 = cos(pkin(11));
t178 = t164 * t166;
t159 = t164 ^ 2;
t161 = t166 ^ 2;
t177 = t161 - t159;
t176 = pkin(9) * m(5) + mrSges(5,3);
t175 = -pkin(10) * m(7) - mrSges(7,3);
t153 = t176 * pkin(3) + Ifges(4,4);
t174 = t153 * t178;
t155 = mrSges(4,2) - t176;
t157 = m(5) * pkin(3) + mrSges(4,1);
t172 = -t164 * t155 + t166 * t157;
t165 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (0.2e1 * t174 + t159 * t152 + t161 * t156 + Ifges(3,2) + ((2 * mrSges(3,3) + t180) * pkin(8))) * t165 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t180) * t165; t181 * t177 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t174; t177 * t153 + t181 * t178 + Ifges(3,4); t166 * Ifges(4,5) - t164 * Ifges(4,6) + Ifges(3,5); t164 * Ifges(4,5) + t166 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + ((t168 + t170) * m(5)) + 0.2e1 * t172 * pkin(2) + t171; mrSges(3,1) + t172; t166 * t155 + t164 * t157 + mrSges(3,2); mrSges(4,3); m(4) + m(5); t169 * m(7) + Ifges(5,1) + Ifges(6,2) - t173 - 2 * t179; Ifges(5,4) + Ifges(6,6); t175 * pkin(5) - Ifges(6,4) + Ifges(5,5); Ifges(5,6) - Ifges(6,5); m(7) * t167 + Ifges(6,1) + Ifges(7,2) + Ifges(5,3) + t162; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) + t175; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
