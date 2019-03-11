% Return the minimum parameter vector for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRR13_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR13_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t202 = sin(pkin(7));
t219 = m(4) * pkin(10) + mrSges(4,3);
t228 = t202 * t219;
t205 = cos(pkin(7));
t227 = t205 * t219;
t206 = (m(6) + m(7));
t226 = (pkin(10) * mrSges(4,3));
t207 = m(3) + m(4);
t225 = t207 * pkin(9);
t212 = (pkin(5) ^ 2);
t224 = (t212 * m(7) + Ifges(6,2));
t223 = 2 * pkin(12) * mrSges(7,3) + Ifges(7,2);
t201 = sin(pkin(13));
t204 = cos(pkin(13));
t222 = t201 * t204;
t213 = (pkin(4) ^ 2);
t191 = (t213 * t206 + Ifges(4,2) + Ifges(5,3));
t221 = Ifges(5,4) * t222;
t220 = 2 * pkin(11) * mrSges(6,3) + t224;
t218 = pkin(12) * m(7) + mrSges(7,3);
t217 = t191 + 2 * t226;
t216 = -pkin(11) * t206 - mrSges(6,3);
t211 = (pkin(10) ^ 2);
t215 = t211 * m(4) + t217;
t214 = pkin(2) ^ 2;
t210 = pkin(11) ^ 2;
t209 = pkin(12) ^ 2;
t203 = sin(pkin(6));
t198 = t205 ^ 2;
t197 = t204 ^ 2;
t194 = t201 ^ 2;
t192 = t211 * t198 + t214;
t190 = pkin(4) * t216 + Ifges(5,5);
t189 = mrSges(3,3) + t227;
t188 = t210 * t206 + Ifges(5,1) + t220;
t187 = Ifges(5,2) + (t210 + t213) * t206 + t220;
t1 = [pkin(1) ^ 2 * t207 + Ifges(2,3) + (t192 * m(4) + Ifges(3,2) + t217 * t198 + (0.2e1 * t189 + t225) * pkin(9)) * t203 ^ 2; pkin(1) * t207 + mrSges(2,1); mrSges(2,2) + (-t189 - t225) * t203; -t198 * t191 + Ifges(3,1) - Ifges(3,2) + (-0.2e1 * t198 + 0.2e1) * t226 + (-t192 + t211) * m(4) + t191; pkin(2) * t228 + Ifges(3,4); -pkin(2) * t227 + Ifges(3,5); t202 * t205 * t215 + Ifges(3,6); t202 ^ 2 * t215 + t214 * m(4) + Ifges(3,3); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t228; t194 * t187 + t197 * t188 + Ifges(4,1) - t191 - 0.2e1 * t221; t201 * Ifges(5,6) - t204 * t190 + Ifges(4,4); Ifges(4,5) + (t197 - t194) * Ifges(5,4) + (-t187 + t188) * t222; -t204 * Ifges(5,6) - t201 * t190 + Ifges(4,6); t197 * t187 + t194 * t188 + Ifges(4,3) + 0.2e1 * t221; mrSges(4,1); mrSges(4,2); pkin(4) * t206 + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t216; m(5) + t206; m(7) * t209 + Ifges(6,1) + t223 - t224; pkin(5) * t218 + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t209 + t212) * m(7) + t223; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t218; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
