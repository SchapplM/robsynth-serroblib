% Return the minimum parameter vector for
% S6RRRRRR9
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRR9_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t199 = sin(pkin(7));
t202 = (m(6) + m(7));
t196 = (m(5) + t202);
t190 = (m(4) + t196);
t215 = pkin(10) * t190 + mrSges(4,3);
t224 = t199 * t215;
t201 = cos(pkin(7));
t223 = t201 * t215;
t222 = (pkin(10) * mrSges(4,3));
t188 = m(3) + t190;
t221 = t188 * pkin(9);
t209 = (pkin(3) ^ 2);
t185 = (t209 * t196 + Ifges(4,2));
t208 = (pkin(4) ^ 2);
t220 = (t208 * t202 + Ifges(5,2));
t203 = (pkin(13) ^ 2);
t207 = (pkin(5) ^ 2);
t219 = (Ifges(6,2) + (t203 + t207) * m(7));
t218 = 2 * pkin(11) * mrSges(5,3) + t220;
t217 = -pkin(13) * m(7) - mrSges(7,3);
t216 = t185 + 2 * t222;
t214 = pkin(11) * t196 + mrSges(5,3);
t187 = (mrSges(6,3) - t217);
t213 = pkin(12) * t202 + t187;
t212 = 2 * pkin(13) * mrSges(7,3) + 2 * pkin(12) * t187 + Ifges(7,2) + t219;
t206 = (pkin(10) ^ 2);
t211 = t206 * t190 + t216;
t210 = pkin(2) ^ 2;
t205 = pkin(11) ^ 2;
t204 = pkin(12) ^ 2;
t200 = sin(pkin(6));
t195 = t201 ^ 2;
t184 = t206 * t195 + t210;
t182 = mrSges(3,3) + t223;
t1 = [pkin(1) ^ 2 * t188 + Ifges(2,3) + (t184 * t190 + Ifges(3,2) + t216 * t195 + (0.2e1 * t182 + t221) * pkin(9)) * t200 ^ 2; pkin(1) * t188 + mrSges(2,1); mrSges(2,2) + (-t182 - t221) * t200; -t195 * t185 + Ifges(3,1) - Ifges(3,2) + (-t184 + t206) * t190 + (-0.2e1 * t195 + 0.2e1) * t222 + t185; pkin(2) * t224 + Ifges(3,4); -pkin(2) * t223 + Ifges(3,5); t199 * t201 * t211 + Ifges(3,6); t199 ^ 2 * t211 + t210 * t190 + Ifges(3,3); pkin(2) * t190 + mrSges(3,1); mrSges(3,2) - t224; t205 * t196 + Ifges(4,1) - t185 + t218; pkin(3) * t214 + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t205 + t209) * t196 + t218; pkin(3) * t196 + mrSges(4,1); mrSges(4,2) - t214; t204 * t202 + Ifges(5,1) + t212 - t220; pkin(4) * t213 + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t204 + t208) * t202 + t212; pkin(4) * t202 + mrSges(5,1); mrSges(5,2) - t213; m(7) * t203 + Ifges(6,1) - t219; Ifges(6,4); pkin(5) * t217 + Ifges(6,5); Ifges(6,6); t207 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
