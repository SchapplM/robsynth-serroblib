% Return the minimum parameter vector for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRP6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t169 = (m(5) + m(6));
t172 = (pkin(9) ^ 2);
t173 = (pkin(4) ^ 2);
t183 = (t173 * m(6) + Ifges(5,2));
t179 = 2 * pkin(9) * mrSges(5,3) + t183;
t156 = t172 * t169 + Ifges(4,1) + t179;
t174 = (pkin(3) ^ 2);
t160 = t174 * t169 + Ifges(4,2);
t186 = t156 - t160;
t185 = (m(3) * pkin(8));
t184 = (-Ifges(6,2) - Ifges(7,3));
t166 = sin(pkin(11));
t168 = cos(pkin(11));
t182 = t166 * t168;
t161 = t166 ^ 2;
t163 = t168 ^ 2;
t181 = t163 - t161;
t180 = 2 * pkin(10) * mrSges(6,3) - t184;
t178 = pkin(10) * m(6) + mrSges(6,3);
t176 = pkin(9) * t169 + mrSges(5,3);
t157 = t176 * pkin(3) + Ifges(4,4);
t177 = t157 * t182;
t158 = mrSges(4,2) - t176;
t159 = pkin(3) * t169 + mrSges(4,1);
t175 = -t166 * t158 + t168 * t159;
t171 = pkin(10) ^ 2;
t167 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (0.2e1 * t177 + t161 * t156 + t163 * t160 + Ifges(3,2) + ((2 * mrSges(3,3) + t185) * pkin(8))) * t167 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t185) * t167; t186 * t181 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t177; t181 * t157 + t186 * t182 + Ifges(3,4); t168 * Ifges(4,5) - t166 * Ifges(4,6) + Ifges(3,5); t166 * Ifges(4,5) + t168 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + ((t172 + t174) * t169) + 0.2e1 * t175 * pkin(2) + t179; mrSges(3,1) + t175; t168 * t158 + t166 * t159 + mrSges(3,2); mrSges(4,3); m(4) + t169; t171 * m(6) + Ifges(5,1) + t180 - t183; t178 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t171 + t173) * m(6) + t180; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t178; Ifges(6,1) + Ifges(7,1) + t184; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
