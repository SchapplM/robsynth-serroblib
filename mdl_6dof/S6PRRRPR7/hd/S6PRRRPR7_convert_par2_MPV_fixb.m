% Return the minimum parameter vector for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% MPV [29x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRRPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t179 = m(4) + m(5);
t192 = (t179 * pkin(9));
t191 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t176 = sin(pkin(13));
t178 = cos(pkin(13));
t190 = t176 * t178;
t183 = (pkin(5) ^ 2);
t189 = (t183 * m(7) + Ifges(5,2) + Ifges(6,3));
t188 = Ifges(6,4) * t190;
t187 = pkin(10) * m(5) + mrSges(5,3);
t186 = -pkin(11) * m(7) - mrSges(7,3);
t185 = 2 * pkin(10) * mrSges(5,3) + t189;
t184 = (pkin(3) ^ 2);
t182 = pkin(10) ^ 2;
t181 = pkin(11) ^ 2;
t177 = sin(pkin(7));
t173 = t178 ^ 2;
t171 = t176 ^ 2;
t170 = t186 * pkin(5) + Ifges(6,5);
t169 = m(7) * t181 + Ifges(6,1) + t191;
t168 = Ifges(6,2) + (t181 + t183) * m(7) + t191;
t1 = [m(2) + m(3) + t179; pkin(2) ^ 2 * t179 + Ifges(3,3) + (t184 * m(5) + Ifges(4,2) + (2 * mrSges(4,3) + t192) * pkin(9)) * t177 ^ 2; pkin(2) * t179 + mrSges(3,1); mrSges(3,2) + (-mrSges(4,3) - t192) * t177; Ifges(4,1) - Ifges(4,2) + (t182 - t184) * m(5) + t185; t187 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t182 + t184) * m(5) + t185; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t187; t171 * t168 + t173 * t169 + Ifges(5,1) - 0.2e1 * t188 - t189; t176 * Ifges(6,6) - t178 * t170 + Ifges(5,4); Ifges(5,5) + (t173 - t171) * Ifges(6,4) + (-t168 + t169) * t190; -t178 * Ifges(6,6) - t176 * t170 + Ifges(5,6); t173 * t168 + t171 * t169 + Ifges(5,3) + 0.2e1 * t188; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t186; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
