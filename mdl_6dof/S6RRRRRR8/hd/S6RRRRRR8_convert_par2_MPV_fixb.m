% Return the minimum parameter vector for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRR8_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t202 = sin(pkin(7));
t205 = (m(6) + m(7));
t199 = (m(5) + t205);
t194 = (m(4) + t199);
t219 = pkin(10) * t194 + mrSges(4,3);
t228 = t202 * t219;
t204 = cos(pkin(7));
t227 = t204 * t219;
t226 = (pkin(10) * mrSges(4,3));
t192 = m(3) + t194;
t225 = t192 * pkin(9);
t213 = (pkin(3) ^ 2);
t191 = (t213 * t199 + Ifges(4,2));
t208 = (pkin(12) ^ 2);
t212 = (pkin(4) ^ 2);
t224 = (Ifges(5,2) + (t208 + t212) * t205);
t211 = (pkin(5) ^ 2);
t223 = (t211 * m(7) + Ifges(6,2));
t222 = 2 * pkin(13) * mrSges(7,3) + Ifges(7,2);
t221 = pkin(13) * m(7) + mrSges(7,3);
t220 = t191 + 2 * t226;
t218 = -pkin(12) * t205 - mrSges(6,3);
t189 = (mrSges(5,3) - t218);
t217 = pkin(11) * t199 + t189;
t216 = 2 * pkin(12) * mrSges(6,3) + 2 * pkin(11) * t189 + t223 + t224;
t210 = (pkin(10) ^ 2);
t215 = t210 * t194 + t220;
t214 = pkin(2) ^ 2;
t209 = pkin(11) ^ 2;
t207 = pkin(13) ^ 2;
t203 = sin(pkin(6));
t198 = t204 ^ 2;
t190 = t210 * t198 + t214;
t186 = mrSges(3,3) + t227;
t1 = [pkin(1) ^ 2 * t192 + Ifges(2,3) + (t190 * t194 + Ifges(3,2) + t220 * t198 + (0.2e1 * t186 + t225) * pkin(9)) * t203 ^ 2; pkin(1) * t192 + mrSges(2,1); mrSges(2,2) + (-t186 - t225) * t203; -t198 * t191 + Ifges(3,1) - Ifges(3,2) + (-t190 + t210) * t194 + (-0.2e1 * t198 + 0.2e1) * t226 + t191; pkin(2) * t228 + Ifges(3,4); -pkin(2) * t227 + Ifges(3,5); t202 * t204 * t215 + Ifges(3,6); t202 ^ 2 * t215 + t214 * t194 + Ifges(3,3); pkin(2) * t194 + mrSges(3,1); mrSges(3,2) - t228; t209 * t199 + Ifges(4,1) - t191 + t216; pkin(3) * t217 + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t209 + t213) * t199 + t216; pkin(3) * t199 + mrSges(4,1); mrSges(4,2) - t217; t208 * t205 + Ifges(5,1) - t224; Ifges(5,4); pkin(4) * t218 + Ifges(5,5); Ifges(5,6); t212 * t205 + Ifges(5,3); pkin(4) * t205 + mrSges(5,1); mrSges(5,2); m(7) * t207 + Ifges(6,1) + t222 - t223; pkin(5) * t221 + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t207 + t211) * m(7) + t222; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t221; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
