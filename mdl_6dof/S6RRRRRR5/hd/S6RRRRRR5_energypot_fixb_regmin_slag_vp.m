% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR5
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:04:46
% EndTime: 2019-03-10 04:04:46
% DurationCPUTime: 0.13s
% Computational Cost: add. (136->52), mult. (194->95), div. (0->0), fcn. (248->14), ass. (0->33)
t271 = sin(pkin(6));
t275 = sin(qJ(2));
t289 = t271 * t275;
t276 = sin(qJ(1));
t288 = t271 * t276;
t278 = cos(qJ(3));
t287 = t271 * t278;
t279 = cos(qJ(2));
t286 = t271 * t279;
t280 = cos(qJ(1));
t285 = t271 * t280;
t284 = t276 * t275;
t283 = t276 * t279;
t282 = t280 * t275;
t281 = t280 * t279;
t270 = qJ(3) + qJ(4);
t277 = cos(qJ(6));
t274 = sin(qJ(3));
t273 = sin(qJ(6));
t272 = cos(pkin(6));
t269 = qJ(5) + t270;
t268 = cos(t270);
t267 = sin(t270);
t266 = cos(t269);
t265 = sin(t269);
t264 = -t272 * t284 + t281;
t263 = t272 * t283 + t282;
t262 = t272 * t282 + t283;
t261 = -t272 * t281 + t284;
t260 = t272 * t265 + t266 * t289;
t259 = t264 * t266 + t265 * t288;
t258 = t262 * t266 - t265 * t285;
t1 = [0, -g(1) * t280 - g(2) * t276, g(1) * t276 - g(2) * t280, 0, 0, 0, 0, 0, -g(1) * t264 - g(2) * t262 - g(3) * t289, g(1) * t263 + g(2) * t261 - g(3) * t286, 0, 0, 0, 0, 0, -g(1) * (t264 * t278 + t274 * t288) - g(2) * (t262 * t278 - t274 * t285) - g(3) * (t272 * t274 + t275 * t287) -g(1) * (-t264 * t274 + t276 * t287) - g(2) * (-t262 * t274 - t278 * t285) - g(3) * (t272 * t278 - t274 * t289) 0, 0, 0, 0, 0, -g(1) * (t264 * t268 + t267 * t288) - g(2) * (t262 * t268 - t267 * t285) - g(3) * (t272 * t267 + t268 * t289) -g(1) * (-t264 * t267 + t268 * t288) - g(2) * (-t262 * t267 - t268 * t285) - g(3) * (-t267 * t289 + t272 * t268) 0, 0, 0, 0, 0, -g(1) * t259 - g(2) * t258 - g(3) * t260, -g(1) * (-t264 * t265 + t266 * t288) - g(2) * (-t262 * t265 - t266 * t285) - g(3) * (-t265 * t289 + t272 * t266) 0, 0, 0, 0, 0, -g(1) * (t259 * t277 + t263 * t273) - g(2) * (t258 * t277 + t261 * t273) - g(3) * (t260 * t277 - t273 * t286) -g(1) * (-t259 * t273 + t263 * t277) - g(2) * (-t258 * t273 + t261 * t277) - g(3) * (-t260 * t273 - t277 * t286);];
U_reg  = t1;
