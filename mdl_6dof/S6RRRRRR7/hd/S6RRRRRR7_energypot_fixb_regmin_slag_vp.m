% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x38]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:41:05
% EndTime: 2019-03-10 04:41:05
% DurationCPUTime: 0.12s
% Computational Cost: add. (132->52), mult. (246->95), div. (0->0), fcn. (320->14), ass. (0->32)
t265 = sin(pkin(6));
t269 = sin(qJ(2));
t282 = t265 * t269;
t272 = cos(qJ(3));
t281 = t265 * t272;
t273 = cos(qJ(2));
t280 = t265 * t273;
t274 = cos(qJ(1));
t279 = t265 * t274;
t270 = sin(qJ(1));
t278 = t270 * t269;
t277 = t270 * t273;
t276 = t274 * t269;
t275 = t274 * t273;
t264 = qJ(4) + qJ(5);
t271 = cos(qJ(4));
t268 = sin(qJ(3));
t267 = sin(qJ(4));
t266 = cos(pkin(6));
t263 = qJ(6) + t264;
t262 = cos(t264);
t261 = sin(t264);
t260 = cos(t263);
t259 = sin(t263);
t258 = -t266 * t278 + t275;
t257 = t266 * t277 + t276;
t256 = t266 * t276 + t277;
t255 = -t266 * t275 + t278;
t254 = t266 * t268 + t269 * t281;
t253 = t270 * t265 * t268 + t258 * t272;
t252 = t256 * t272 - t268 * t279;
t1 = [0, -g(1) * t274 - g(2) * t270, g(1) * t270 - g(2) * t274, 0, 0, 0, 0, 0, -g(1) * t258 - g(2) * t256 - g(3) * t282, g(1) * t257 + g(2) * t255 - g(3) * t280, 0, 0, 0, 0, 0, -g(1) * t253 - g(2) * t252 - g(3) * t254, -g(1) * (-t258 * t268 + t270 * t281) - g(2) * (-t256 * t268 - t272 * t279) - g(3) * (t266 * t272 - t268 * t282) 0, 0, 0, 0, 0, -g(1) * (t253 * t271 + t257 * t267) - g(2) * (t252 * t271 + t255 * t267) - g(3) * (t254 * t271 - t267 * t280) -g(1) * (-t253 * t267 + t257 * t271) - g(2) * (-t252 * t267 + t255 * t271) - g(3) * (-t254 * t267 - t271 * t280) 0, 0, 0, 0, 0, -g(1) * (t253 * t262 + t257 * t261) - g(2) * (t252 * t262 + t255 * t261) - g(3) * (t254 * t262 - t261 * t280) -g(1) * (-t253 * t261 + t257 * t262) - g(2) * (-t252 * t261 + t255 * t262) - g(3) * (-t254 * t261 - t262 * t280) 0, 0, 0, 0, 0, -g(1) * (t253 * t260 + t257 * t259) - g(2) * (t252 * t260 + t255 * t259) - g(3) * (t254 * t260 - t259 * t280) -g(1) * (-t253 * t259 + t257 * t260) - g(2) * (-t252 * t259 + t255 * t260) - g(3) * (-t254 * t259 - t260 * t280);];
U_reg  = t1;
