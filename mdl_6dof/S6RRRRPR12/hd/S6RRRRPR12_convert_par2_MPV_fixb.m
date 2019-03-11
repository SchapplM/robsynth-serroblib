% Return the minimum parameter vector for
% S6RRRRPR12
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
% MPV [33x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRPR12_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR12_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t207 = sin(pkin(7));
t211 = (m(4) + m(5));
t222 = pkin(10) * t211 + mrSges(4,3);
t233 = t222 * t207;
t210 = cos(pkin(7));
t232 = t222 * t210;
t231 = (pkin(10) * mrSges(4,3));
t203 = m(3) + t211;
t230 = t203 * pkin(9);
t217 = (pkin(3) ^ 2);
t197 = (t217 * m(5) + Ifges(4,2));
t229 = 2 * pkin(12) * mrSges(7,3) + Ifges(7,2);
t206 = sin(pkin(13));
t209 = cos(pkin(13));
t228 = t206 * t209;
t213 = pkin(12) ^ 2;
t190 = m(7) * t213 + Ifges(6,1) + t229;
t216 = (pkin(5) ^ 2);
t195 = t216 * m(7) + Ifges(6,2);
t198 = t206 ^ 2;
t201 = t209 ^ 2;
t227 = t198 * t190 + t201 * t195 + Ifges(5,2);
t226 = pkin(11) * m(5) + mrSges(5,3);
t225 = pkin(12) * m(7) + mrSges(7,3);
t192 = t225 * pkin(5) + Ifges(6,4);
t224 = t192 * t228;
t223 = t197 + 2 * t231;
t221 = (2 * pkin(11) * mrSges(5,3)) + 0.2e1 * t224 + t227;
t194 = mrSges(6,2) - t225;
t196 = m(7) * pkin(5) + mrSges(6,1);
t220 = -t206 * t194 + t209 * t196;
t215 = (pkin(10) ^ 2);
t219 = t215 * t211 + t223;
t218 = pkin(2) ^ 2;
t214 = pkin(11) ^ 2;
t208 = sin(pkin(6));
t202 = t210 ^ 2;
t193 = t215 * t202 + t218;
t188 = mrSges(3,3) + t232;
t1 = [pkin(1) ^ 2 * t203 + Ifges(2,3) + (t193 * t211 + Ifges(3,2) + t223 * t202 + (0.2e1 * t188 + t230) * pkin(9)) * t208 ^ 2; pkin(1) * t203 + mrSges(2,1); mrSges(2,2) + (-t188 - t230) * t208; -t202 * t197 + Ifges(3,1) - Ifges(3,2) + (-t193 + t215) * t211 + (-0.2e1 * t202 + 0.2e1) * t231 + t197; pkin(2) * t233 + Ifges(3,4); -pkin(2) * t232 + Ifges(3,5); t219 * t210 * t207 + Ifges(3,6); t219 * t207 ^ 2 + t218 * t211 + Ifges(3,3); pkin(2) * t211 + mrSges(3,1); mrSges(3,2) - t233; (t214 * m(5)) + Ifges(4,1) - t197 + t221; t226 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + ((t214 + t217) * m(5)) + t221; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t226; t201 * t190 + t198 * t195 + Ifges(5,1) - 0.4e1 * t224 - t227; Ifges(5,4) + (t201 - t198) * t192 + (t190 - t195) * t228; t209 * Ifges(6,5) - t206 * Ifges(6,6) + Ifges(5,5); t206 * Ifges(6,5) + t209 * Ifges(6,6) + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t213 + t216) * m(7)) + 0.2e1 * t220 * pkin(4) + t229; mrSges(5,1) + t220; t209 * t194 + t206 * t196 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
