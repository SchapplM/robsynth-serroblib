% Return the minimum parameter vector for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRR14_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_convert_par2_MPV_fixb: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t259 = sin(pkin(8));
t263 = cos(pkin(8));
t265 = (m(6) + m(7));
t271 = (pkin(4) ^ 2);
t246 = (t271 * t265 + Ifges(5,2));
t255 = m(5) + t265;
t298 = (pkin(11) * t255);
t274 = t246 + (2 * mrSges(5,3) + t298) * pkin(11);
t237 = t274 * t263 * t259 + Ifges(4,6);
t280 = mrSges(5,3) + t298;
t296 = t280 * t263;
t241 = -pkin(3) * t296 + Ifges(4,5);
t258 = sin(pkin(14));
t262 = cos(pkin(14));
t275 = t237 * t262 + t241 * t258;
t260 = sin(pkin(7));
t264 = cos(pkin(7));
t287 = t260 * t264;
t299 = t275 * t287;
t297 = t280 * t259;
t286 = pkin(3) ^ 2 * t255;
t238 = t274 * t259 ^ 2 + Ifges(4,3) + t286;
t250 = t260 ^ 2;
t254 = t264 ^ 2;
t239 = t274 * t263 ^ 2 + Ifges(4,2) + t286;
t243 = Ifges(4,1) + t274;
t248 = t258 ^ 2;
t252 = t262 ^ 2;
t242 = pkin(3) * t297 + Ifges(4,4);
t281 = t262 * t258 * t242;
t273 = t239 * t252 + t243 * t248 + 0.2e1 * t281;
t295 = t250 * t238 + t273 * t254 + Ifges(3,2);
t294 = (m(3) * pkin(10));
t270 = (pkin(5) ^ 2);
t293 = (t270 * m(7) + Ifges(6,2));
t292 = 2 * pkin(13) * mrSges(7,3) + Ifges(7,2);
t289 = t241 * t262;
t288 = t242 * (t252 - t248);
t285 = 0.2e1 * t299;
t283 = 2 * pkin(12) * mrSges(6,3) + t293;
t282 = m(7) * pkin(13) + mrSges(7,3);
t279 = pkin(12) * t265 + mrSges(6,3);
t277 = (-t239 + t243) * t262;
t268 = pkin(12) ^ 2;
t267 = pkin(13) ^ 2;
t261 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (-0.2e1 * t299 + ((2 * mrSges(3,3) + t294) * pkin(10)) + t295) * t261 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t294) * t261; t248 * t239 + t252 * t243 + Ifges(3,1) - 0.2e1 * t281 + t285 - t295; -t260 * t289 + t264 * t288 + Ifges(3,4) + (t237 * t260 + t264 * t277) * t258; t264 * t289 + t260 * t288 + Ifges(3,5) + (-t237 * t264 + t260 * t277) * t258; Ifges(3,6) + t275 * (t254 - t250) + (-t238 + t273) * t287; t254 * t238 + t273 * t250 + Ifges(3,3) + t285; mrSges(3,1); mrSges(3,2); pkin(3) * t255 + mrSges(4,1); mrSges(4,2) - t297; mrSges(4,3) + t296; m(4) + t255; t268 * t265 + Ifges(5,1) - t246 + t283; t279 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t268 + t271) * t265 + t283; pkin(4) * t265 + mrSges(5,1); mrSges(5,2) - t279; m(7) * t267 + Ifges(6,1) + t292 - t293; t282 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t267 + t270) * m(7) + t292; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t282; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
