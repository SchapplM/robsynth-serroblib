% Return the minimum parameter vector for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRR15_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR15_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t193 = sin(pkin(7));
t210 = m(4) * pkin(10) + mrSges(4,3);
t217 = t210 * t193;
t195 = cos(pkin(7));
t216 = t210 * t195;
t196 = (m(6) + m(7));
t215 = (pkin(10) * mrSges(4,3));
t214 = (pkin(11) * mrSges(6,3));
t197 = m(3) + m(4);
t213 = t197 * pkin(9);
t202 = (pkin(5) ^ 2);
t212 = (t202 * m(7) + Ifges(6,2));
t211 = 2 * pkin(12) * mrSges(7,3) + Ifges(7,2);
t209 = pkin(12) * m(7) + mrSges(7,3);
t192 = 2 * t214;
t200 = (pkin(11) ^ 2);
t203 = (pkin(4) ^ 2);
t206 = (Ifges(4,2) + Ifges(5,3) + (t200 + t203) * t196 + t212);
t183 = t192 + t206;
t208 = t183 + 2 * t215;
t207 = -pkin(11) * t196 - mrSges(6,3);
t201 = (pkin(10) ^ 2);
t205 = t201 * m(4) + t208;
t204 = pkin(2) ^ 2;
t199 = pkin(12) ^ 2;
t194 = sin(pkin(6));
t190 = t195 ^ 2;
t186 = t201 * t190 + t204;
t184 = mrSges(3,3) + t216;
t1 = [pkin(1) ^ 2 * t197 + Ifges(2,3) + (t186 * m(4) + Ifges(3,2) + t208 * t190 + (0.2e1 * t184 + t213) * pkin(9)) * t194 ^ 2; pkin(1) * t197 + mrSges(2,1); mrSges(2,2) + (-t184 - t213) * t194; -t190 * t183 + Ifges(3,1) - Ifges(3,2) + (-0.2e1 * t190 + 0.2e1) * t215 + (-t186 + t201) * m(4) + t183; pkin(2) * t217 + Ifges(3,4); -pkin(2) * t216 + Ifges(3,5); t205 * t195 * t193 + Ifges(3,6); t205 * t193 ^ 2 + t204 * m(4) + Ifges(3,3); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t217; t203 * t196 + Ifges(4,1) + Ifges(5,2) - t206 - 2 * t214; Ifges(4,4) + Ifges(5,6); t207 * pkin(4) - Ifges(5,4) + Ifges(4,5); Ifges(4,6) - Ifges(5,5); t200 * t196 + Ifges(5,1) + Ifges(4,3) + t192 + t212; mrSges(4,1); mrSges(4,2); pkin(4) * t196 + mrSges(5,1); mrSges(5,2) + t207; mrSges(5,3); m(5) + t196; m(7) * t199 + Ifges(6,1) + t211 - t212; t209 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t199 + t202) * m(7) + t211; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t209; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
