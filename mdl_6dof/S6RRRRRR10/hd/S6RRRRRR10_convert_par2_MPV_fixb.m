% Return the minimum parameter vector for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_convert_par2_MPV_fixb: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t239 = sin(pkin(7));
t241 = cos(pkin(8));
t243 = (m(6) + m(7));
t235 = m(5) + t243;
t258 = pkin(12) * t235 + mrSges(5,3);
t270 = t258 * t241;
t220 = mrSges(4,3) + t270;
t228 = m(4) + t235;
t255 = pkin(11) * t228 + t220;
t273 = t255 * t239;
t242 = cos(pkin(7));
t272 = t255 * t242;
t238 = sin(pkin(8));
t271 = t258 * t238;
t250 = (pkin(4) ^ 2);
t226 = t250 * t243 + Ifges(5,2);
t247 = pkin(12) ^ 2;
t269 = (t247 * t235 + t226);
t268 = (pkin(12) * mrSges(5,3));
t267 = pkin(11) * t220;
t227 = m(3) + t228;
t266 = t227 * pkin(10);
t249 = (pkin(5) ^ 2);
t265 = (t249 * m(7) + Ifges(6,2));
t264 = 2 * pkin(14) * mrSges(7,3) + Ifges(7,2);
t262 = 2 * t268;
t233 = t241 ^ 2;
t251 = pkin(3) ^ 2;
t261 = Ifges(4,2) + (t247 * t233 + t251) * t235 + t233 * t226;
t260 = 2 * pkin(13) * mrSges(6,3) + t265;
t259 = pkin(14) * m(7) + mrSges(7,3);
t257 = pkin(13) * t243 + mrSges(6,3);
t219 = t233 * t262 + t261;
t256 = t219 + 0.2e1 * t267;
t254 = t262 + t269;
t248 = pkin(11) ^ 2;
t253 = t248 * t228 + t256;
t252 = pkin(2) ^ 2;
t246 = pkin(13) ^ 2;
t245 = pkin(14) ^ 2;
t240 = sin(pkin(6));
t234 = t242 ^ 2;
t224 = t248 * t234 + t252;
t218 = mrSges(3,3) + t272;
t1 = [pkin(1) ^ 2 * t227 + Ifges(2,3) + (t224 * t228 + Ifges(3,2) + t256 * t234 + (0.2e1 * t218 + t266) * pkin(10)) * t240 ^ 2; pkin(1) * t227 + mrSges(2,1); mrSges(2,2) + (-t218 - t266) * t240; -t234 * t219 + Ifges(3,1) - Ifges(3,2) + (-t224 + t248) * t228 + (-0.2e1 * t234 + 0.2e1) * t267 + t219; pkin(2) * t273 + Ifges(3,4); -pkin(2) * t272 + Ifges(3,5); t253 * t242 * t239 + Ifges(3,6); t253 * t239 ^ 2 + t252 * t228 + Ifges(3,3); pkin(2) * t228 + mrSges(3,1); mrSges(3,2) - t273; Ifges(4,1) + (-0.2e1 * t233 + 0.2e1) * t268 - t261 + t269; pkin(3) * t271 + Ifges(4,4); -pkin(3) * t270 + Ifges(4,5); t254 * t241 * t238 + Ifges(4,6); t254 * t238 ^ 2 + t251 * t235 + Ifges(4,3); pkin(3) * t235 + mrSges(4,1); mrSges(4,2) - t271; (t246 * t243) + Ifges(5,1) - t226 + t260; t257 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t246 + t250) * t243 + t260; pkin(4) * t243 + mrSges(5,1); mrSges(5,2) - t257; m(7) * t245 + Ifges(6,1) + t264 - t265; t259 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t245 + t249) * m(7) + t264; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t259; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
