% Return the minimum parameter vector for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRPR14_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR14_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t197 = sin(pkin(7));
t201 = (m(4) + m(5));
t212 = pkin(10) * t201 + mrSges(4,3);
t223 = t197 * t212;
t200 = cos(pkin(7));
t222 = t200 * t212;
t221 = (pkin(10) * mrSges(4,3));
t193 = m(3) + t201;
t220 = t193 * pkin(9);
t208 = (pkin(3) ^ 2);
t187 = (t208 * m(5) + Ifges(4,2));
t219 = 2 * pkin(12) * mrSges(7,3) + Ifges(7,2);
t196 = sin(pkin(13));
t199 = cos(pkin(13));
t218 = t196 * t199;
t207 = (pkin(5) ^ 2);
t217 = (t207 * m(7) + Ifges(5,2) + Ifges(6,3));
t216 = Ifges(6,4) * t218;
t215 = pkin(11) * m(5) + mrSges(5,3);
t214 = -pkin(12) * m(7) - mrSges(7,3);
t213 = t187 + 2 * t221;
t211 = 2 * pkin(11) * mrSges(5,3) + t217;
t206 = (pkin(10) ^ 2);
t210 = t206 * t201 + t213;
t209 = pkin(2) ^ 2;
t205 = pkin(11) ^ 2;
t204 = pkin(12) ^ 2;
t198 = sin(pkin(6));
t192 = t200 ^ 2;
t191 = t199 ^ 2;
t188 = t196 ^ 2;
t186 = t206 * t192 + t209;
t185 = pkin(5) * t214 + Ifges(6,5);
t184 = m(7) * t204 + Ifges(6,1) + t219;
t183 = Ifges(6,2) + (t204 + t207) * m(7) + t219;
t182 = mrSges(3,3) + t222;
t1 = [pkin(1) ^ 2 * t193 + Ifges(2,3) + (t186 * t201 + Ifges(3,2) + t213 * t192 + (0.2e1 * t182 + t220) * pkin(9)) * t198 ^ 2; pkin(1) * t193 + mrSges(2,1); mrSges(2,2) + (-t182 - t220) * t198; -t192 * t187 + Ifges(3,1) - Ifges(3,2) + (-t186 + t206) * t201 + (-0.2e1 * t192 + 0.2e1) * t221 + t187; pkin(2) * t223 + Ifges(3,4); -pkin(2) * t222 + Ifges(3,5); t197 * t200 * t210 + Ifges(3,6); t197 ^ 2 * t210 + t209 * t201 + Ifges(3,3); pkin(2) * t201 + mrSges(3,1); mrSges(3,2) - t223; t205 * m(5) + Ifges(4,1) - t187 + t211; pkin(3) * t215 + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t205 + t208) * m(5) + t211; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t215; t188 * t183 + t191 * t184 + Ifges(5,1) - 0.2e1 * t216 - t217; t196 * Ifges(6,6) - t199 * t185 + Ifges(5,4); Ifges(5,5) + (t191 - t188) * Ifges(6,4) + (-t183 + t184) * t218; -t199 * Ifges(6,6) - t196 * t185 + Ifges(5,6); t191 * t183 + t188 * t184 + Ifges(5,3) + 0.2e1 * t216; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t214; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
