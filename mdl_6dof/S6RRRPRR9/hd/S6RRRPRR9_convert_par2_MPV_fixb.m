% Return the minimum parameter vector for
% S6RRRPRR9
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
% MPV [33x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRR9_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t203 = sin(pkin(7));
t222 = m(4) * pkin(10) + mrSges(4,3);
t231 = t222 * t203;
t206 = cos(pkin(7));
t230 = t222 * t206;
t207 = (m(6) + m(7));
t229 = (pkin(10) * mrSges(4,3));
t208 = m(3) + m(4);
t228 = t208 * pkin(9);
t213 = (pkin(5) ^ 2);
t227 = (t213 * m(7) + Ifges(6,2));
t226 = 2 * pkin(12) * mrSges(7,3) + Ifges(7,2);
t202 = sin(pkin(13));
t205 = cos(pkin(13));
t225 = t202 * t205;
t211 = (pkin(11) ^ 2);
t223 = 2 * pkin(11) * mrSges(6,3) + t227;
t187 = t211 * t207 + Ifges(5,1) + t223;
t214 = (pkin(4) ^ 2);
t194 = t214 * t207 + Ifges(5,2);
t195 = t202 ^ 2;
t198 = t205 ^ 2;
t224 = t195 * t187 + t198 * t194 + Ifges(4,2);
t221 = pkin(12) * m(7) + mrSges(7,3);
t218 = pkin(11) * t207 + mrSges(6,3);
t190 = t218 * pkin(4) + Ifges(5,4);
t220 = t190 * t225;
t184 = 0.2e1 * t220 + t224;
t219 = t184 + (2 * t229);
t191 = mrSges(5,2) - t218;
t193 = pkin(4) * t207 + mrSges(5,1);
t217 = -t202 * t191 + t205 * t193;
t212 = pkin(10) ^ 2;
t216 = t212 * m(4) + t219;
t215 = pkin(2) ^ 2;
t210 = pkin(12) ^ 2;
t204 = sin(pkin(6));
t199 = t206 ^ 2;
t192 = t212 * t199 + t215;
t188 = mrSges(3,3) + t230;
t1 = [pkin(1) ^ 2 * t208 + Ifges(2,3) + (t192 * m(4) + Ifges(3,2) + t219 * t199 + (0.2e1 * t188 + t228) * pkin(9)) * t204 ^ 2; pkin(1) * t208 + mrSges(2,1); mrSges(2,2) + (-t188 - t228) * t204; -t199 * t184 + Ifges(3,1) - Ifges(3,2) + (-0.2e1 * t199 + 0.2e1) * t229 + (-t192 + t212) * m(4) + t184; pkin(2) * t231 + Ifges(3,4); -pkin(2) * t230 + Ifges(3,5); t216 * t206 * t203 + Ifges(3,6); t216 * t203 ^ 2 + t215 * m(4) + Ifges(3,3); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t231; t198 * t187 + t195 * t194 + Ifges(4,1) - 0.4e1 * t220 - t224; Ifges(4,4) + (t198 - t195) * t190 + (t187 - t194) * t225; t205 * Ifges(5,5) - t202 * Ifges(5,6) + Ifges(4,5); t202 * Ifges(5,5) + t205 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t211 + t214) * t207) + 0.2e1 * t217 * pkin(3) + t223; mrSges(4,1) + t217; t205 * t191 + t202 * t193 + mrSges(4,2); mrSges(5,3); m(5) + t207; m(7) * t210 + Ifges(6,1) + t226 - t227; t221 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t210 + t213) * m(7) + t226; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t221; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
