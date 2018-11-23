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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_convert_par2_MPV_fixb: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t257 = sin(pkin(7));
t259 = cos(pkin(8));
t261 = (m(6) + m(7));
t253 = m(5) + t261;
t276 = pkin(12) * t253 + mrSges(5,3);
t288 = t276 * t259;
t238 = mrSges(4,3) + t288;
t246 = m(4) + t253;
t273 = pkin(11) * t246 + t238;
t291 = t273 * t257;
t260 = cos(pkin(7));
t290 = t273 * t260;
t256 = sin(pkin(8));
t289 = t276 * t256;
t268 = (pkin(4) ^ 2);
t244 = t268 * t261 + Ifges(5,2);
t265 = pkin(12) ^ 2;
t287 = (t253 * t265 + t244);
t286 = (pkin(12) * mrSges(5,3));
t245 = m(3) + t246;
t285 = pkin(10) * t245;
t284 = pkin(11) * t238;
t267 = (pkin(5) ^ 2);
t283 = (t267 * m(7) + Ifges(6,2));
t282 = 2 * pkin(14) * mrSges(7,3) + Ifges(7,2);
t280 = 2 * t286;
t251 = t259 ^ 2;
t269 = pkin(3) ^ 2;
t279 = Ifges(4,2) + (t251 * t265 + t269) * t253 + t251 * t244;
t278 = 2 * pkin(13) * mrSges(6,3) + t283;
t277 = pkin(14) * m(7) + mrSges(7,3);
t275 = pkin(13) * t261 + mrSges(6,3);
t237 = t251 * t280 + t279;
t274 = t237 + 0.2e1 * t284;
t272 = t280 + t287;
t266 = pkin(11) ^ 2;
t271 = t246 * t266 + t274;
t270 = pkin(2) ^ 2;
t264 = pkin(13) ^ 2;
t263 = pkin(14) ^ 2;
t258 = sin(pkin(6));
t252 = t260 ^ 2;
t242 = t252 * t266 + t270;
t236 = mrSges(3,3) + t290;
t1 = [pkin(1) ^ 2 * t245 + Ifges(2,3) + (t242 * t246 + Ifges(3,2) + t274 * t252 + (0.2e1 * t236 + t285) * pkin(10)) * t258 ^ 2; pkin(1) * t245 + mrSges(2,1); mrSges(2,2) + (-t236 - t285) * t258; -t252 * t237 + Ifges(3,1) - Ifges(3,2) + (-t242 + t266) * t246 + (-0.2e1 * t252 + 0.2e1) * t284 + t237; pkin(2) * t291 + Ifges(3,4); -pkin(2) * t290 + Ifges(3,5); t271 * t260 * t257 + Ifges(3,6); t271 * t257 ^ 2 + t270 * t246 + Ifges(3,3); pkin(2) * t246 + mrSges(3,1); mrSges(3,2) - t291; Ifges(4,1) + (-0.2e1 * t251 + 0.2e1) * t286 - t279 + t287; pkin(3) * t289 + Ifges(4,4); -pkin(3) * t288 + Ifges(4,5); t272 * t259 * t256 + Ifges(4,6); t272 * t256 ^ 2 + t269 * t253 + Ifges(4,3); pkin(3) * t253 + mrSges(4,1); mrSges(4,2) - t289; (t261 * t264) + Ifges(5,1) - t244 + t278; t275 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t264 + t268) * t261 + t278; pkin(4) * t261 + mrSges(5,1); mrSges(5,2) - t275; m(7) * t263 + Ifges(6,1) + t282 - t283; t277 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t263 + t267) * m(7) + t282; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t277; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
