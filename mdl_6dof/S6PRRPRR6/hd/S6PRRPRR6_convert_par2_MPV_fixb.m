% Return the minimum parameter vector for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRPRR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t197 = (m(4) * pkin(9));
t183 = (m(6) + m(7));
t187 = (pkin(5) ^ 2);
t196 = (t187 * m(7) + Ifges(6,2));
t195 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t180 = sin(pkin(13));
t182 = cos(pkin(13));
t194 = t180 * t182;
t193 = Ifges(5,4) * t194;
t192 = 2 * pkin(10) * mrSges(6,3) + t196;
t191 = pkin(11) * m(7) + mrSges(7,3);
t190 = -pkin(10) * t183 - mrSges(6,3);
t188 = (pkin(4) ^ 2);
t189 = t188 * t183 + Ifges(4,2) + Ifges(5,3);
t186 = pkin(10) ^ 2;
t185 = pkin(11) ^ 2;
t181 = sin(pkin(7));
t177 = t182 ^ 2;
t175 = t180 ^ 2;
t174 = t190 * pkin(4) + Ifges(5,5);
t173 = t186 * t183 + Ifges(5,1) + t192;
t172 = Ifges(5,2) + (t186 + t188) * t183 + t192;
t1 = [m(2) + m(3) + m(4); pkin(2) ^ 2 * m(4) + Ifges(3,3) + ((2 * mrSges(4,3) + t197) * pkin(9) + t189) * t181 ^ 2; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) + (-mrSges(4,3) - t197) * t181; t175 * t172 + t177 * t173 + Ifges(4,1) - t189 - 0.2e1 * t193; t180 * Ifges(5,6) - t182 * t174 + Ifges(4,4); Ifges(4,5) + (t177 - t175) * Ifges(5,4) + (-t172 + t173) * t194; -t182 * Ifges(5,6) - t180 * t174 + Ifges(4,6); t177 * t172 + t175 * t173 + Ifges(4,3) + 0.2e1 * t193; mrSges(4,1); mrSges(4,2); pkin(4) * t183 + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t190; m(5) + t183; m(7) * t185 + Ifges(6,1) + t195 - t196; t191 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t185 + t187) * m(7) + t195; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t191; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
