% Calculate minimal parameter regressor of potential energy for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:23
% EndTime: 2019-03-08 19:05:23
% DurationCPUTime: 0.14s
% Computational Cost: add. (214->58), mult. (557->111), div. (0->0), fcn. (733->16), ass. (0->43)
t257 = sin(pkin(12));
t263 = cos(pkin(6));
t278 = t257 * t263;
t258 = sin(pkin(7));
t259 = sin(pkin(6));
t277 = t258 * t259;
t276 = t258 * t263;
t261 = cos(pkin(12));
t275 = t259 * t261;
t262 = cos(pkin(7));
t274 = t259 * t262;
t260 = cos(pkin(13));
t273 = t260 * t262;
t272 = t261 * t263;
t256 = sin(pkin(13));
t249 = -t257 * t256 + t260 * t272;
t271 = -t249 * t262 + t258 * t275;
t251 = -t261 * t256 - t260 * t278;
t270 = t251 * t262 + t257 * t277;
t269 = cos(qJ(3));
t268 = cos(qJ(4));
t267 = cos(qJ(5));
t266 = sin(qJ(3));
t265 = sin(qJ(4));
t264 = sin(qJ(5));
t255 = qJ(5) + qJ(6);
t254 = cos(t255);
t253 = sin(t255);
t252 = -t256 * t278 + t261 * t260;
t250 = t256 * t272 + t257 * t260;
t248 = -t260 * t277 + t263 * t262;
t247 = -t251 * t258 + t257 * t274;
t246 = -t249 * t258 - t261 * t274;
t245 = t266 * t276 + (t256 * t269 + t266 * t273) * t259;
t244 = -t269 * t276 + (t256 * t266 - t269 * t273) * t259;
t243 = t245 * t268 + t248 * t265;
t242 = t252 * t269 + t270 * t266;
t241 = t252 * t266 - t270 * t269;
t240 = t250 * t269 - t271 * t266;
t239 = t250 * t266 + t271 * t269;
t238 = t242 * t268 + t247 * t265;
t237 = t240 * t268 + t246 * t265;
t1 = [-g(3) * qJ(1), -g(1) * (t257 * t259 * qJ(2) + t261 * pkin(1)) - g(2) * (t257 * pkin(1) - qJ(2) * t275) - g(3) * (t263 * qJ(2) + qJ(1)) 0, -g(1) * t242 - g(2) * t240 - g(3) * t245, g(1) * t241 + g(2) * t239 + g(3) * t244, 0, 0, 0, 0, 0, -g(1) * t238 - g(2) * t237 - g(3) * t243, -g(1) * (-t242 * t265 + t247 * t268) - g(2) * (-t240 * t265 + t246 * t268) - g(3) * (-t245 * t265 + t248 * t268) 0, 0, 0, 0, 0, -g(1) * (t238 * t267 + t241 * t264) - g(2) * (t237 * t267 + t239 * t264) - g(3) * (t243 * t267 + t244 * t264) -g(1) * (-t238 * t264 + t241 * t267) - g(2) * (-t237 * t264 + t239 * t267) - g(3) * (-t243 * t264 + t244 * t267) 0, 0, 0, 0, 0, -g(1) * (t238 * t254 + t241 * t253) - g(2) * (t237 * t254 + t239 * t253) - g(3) * (t243 * t254 + t244 * t253) -g(1) * (-t238 * t253 + t241 * t254) - g(2) * (-t237 * t253 + t239 * t254) - g(3) * (-t243 * t253 + t244 * t254);];
U_reg  = t1;
