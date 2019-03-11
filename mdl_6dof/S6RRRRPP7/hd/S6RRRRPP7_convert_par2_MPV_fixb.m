% Return the minimum parameter vector for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRPP7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t169 = (m(4) + m(5));
t158 = m(3) + t169;
t184 = (t158 * pkin(8));
t173 = (pkin(3) ^ 2);
t183 = (t173 * m(5) + Ifges(4,2));
t161 = sin(pkin(11));
t163 = cos(pkin(11));
t182 = t161 * t163;
t181 = 2 * pkin(9) * mrSges(4,3) + t183;
t155 = t161 ^ 2;
t157 = t163 ^ 2;
t165 = Ifges(6,2) + Ifges(7,3);
t168 = Ifges(6,1) + Ifges(7,1);
t180 = t155 * t168 + t157 * t165 + Ifges(5,2);
t179 = pkin(10) * m(5) + mrSges(5,3);
t167 = Ifges(6,4) - Ifges(7,5);
t178 = t167 * t182;
t177 = pkin(9) * t169 + mrSges(4,3);
t176 = (2 * pkin(10) * mrSges(5,3)) + 0.2e1 * t178 + t180;
t175 = t163 * mrSges(6,1) - t161 * mrSges(6,2);
t174 = (pkin(2) ^ 2);
t172 = pkin(9) ^ 2;
t171 = pkin(10) ^ 2;
t166 = Ifges(6,5) + Ifges(7,4);
t164 = Ifges(6,6) - Ifges(7,6);
t162 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * t158 + Ifges(2,3) + (t174 * t169 + Ifges(3,2) + (2 * mrSges(3,3) + t184) * pkin(8)) * t162 ^ 2; pkin(1) * t158 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t184) * t162; Ifges(3,1) - Ifges(3,2) + (t172 - t174) * t169 + t181; t177 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t172 + t174) * t169 + t181; pkin(2) * t169 + mrSges(3,1); mrSges(3,2) - t177; (t171 * m(5)) + Ifges(4,1) + t176 - t183; t179 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + ((t171 + t173) * m(5)) + t176; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t179; t155 * t165 + t157 * t168 + Ifges(5,1) - 0.4e1 * t178 - t180; Ifges(5,4) + (t157 - t155) * t167 + (-t165 + t168) * t182; -t161 * t164 + t163 * t166 + Ifges(5,5); t161 * t166 + t163 * t164 + Ifges(5,6); 0.2e1 * pkin(4) * t175 + Ifges(7,2) + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t175; t161 * mrSges(6,1) + t163 * mrSges(6,2) + mrSges(5,2); mrSges(6,3); m(6); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
