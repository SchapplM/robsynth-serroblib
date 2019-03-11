% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:29:33
% EndTime: 2019-03-09 11:29:33
% DurationCPUTime: 0.18s
% Computational Cost: add. (155->72), mult. (338->110), div. (0->0), fcn. (415->12), ass. (0->41)
t260 = sin(qJ(1));
t280 = g(1) * t260;
t263 = cos(qJ(1));
t279 = g(2) * t263;
t255 = sin(pkin(6));
t259 = sin(qJ(2));
t278 = t255 * t259;
t277 = t255 * t260;
t262 = cos(qJ(2));
t276 = t255 * t262;
t275 = t255 * t263;
t274 = t260 * t259;
t273 = t260 * t262;
t272 = t263 * t259;
t271 = t263 * t262;
t257 = cos(pkin(6));
t270 = pkin(2) * t278 + t257 * pkin(8) + pkin(7);
t269 = -t279 + t280;
t240 = -t257 * t271 + t274;
t241 = t257 * t272 + t273;
t268 = t260 * pkin(1) + t241 * pkin(2) + t240 * qJ(3);
t242 = t257 * t273 + t272;
t243 = -t257 * t274 + t271;
t267 = t263 * pkin(1) + t243 * pkin(2) + pkin(8) * t277 + t242 * qJ(3);
t258 = sin(qJ(4));
t261 = cos(qJ(4));
t231 = -t242 * t261 + t258 * t277;
t233 = t240 * t261 + t258 * t275;
t238 = t257 * t258 + t261 * t276;
t266 = g(1) * t231 - g(2) * t233 + g(3) * t238;
t265 = -g(1) * t242 - g(2) * t240 + g(3) * t276;
t264 = g(1) * t243 + g(2) * t241 + g(3) * t278;
t256 = cos(pkin(11));
t254 = sin(pkin(11));
t253 = pkin(11) + qJ(6);
t249 = cos(t253);
t248 = sin(t253);
t239 = t257 * t261 - t258 * t276;
t234 = t240 * t258 - t261 * t275;
t232 = t242 * t258 + t261 * t277;
t1 = [0, -g(1) * t263 - g(2) * t260, t269, 0, 0, 0, 0, 0, -t264, -t265, -g(3) * t257 - t269 * t255, t264, t265, -g(1) * t267 - g(2) * (-pkin(8) * t275 + t268) - g(3) * (-qJ(3) * t276 + t270) 0, 0, 0, 0, 0, -g(1) * t232 - g(2) * t234 - g(3) * t239, t266, -g(1) * (t232 * t256 + t243 * t254) - g(2) * (t234 * t256 + t241 * t254) - g(3) * (t239 * t256 + t254 * t278) -g(1) * (-t232 * t254 + t243 * t256) - g(2) * (-t234 * t254 + t241 * t256) - g(3) * (-t239 * t254 + t256 * t278) -t266, -g(1) * (t232 * pkin(4) + t243 * pkin(9) + t231 * qJ(5) + t267) - g(2) * (t234 * pkin(4) + t241 * pkin(9) - t233 * qJ(5) + t268) - g(3) * (t257 * pkin(3) + t239 * pkin(4) + t238 * qJ(5) + t270) + (-pkin(3) * t280 - g(3) * (pkin(9) * t259 - qJ(3) * t262) - (-pkin(3) - pkin(8)) * t279) * t255, 0, 0, 0, 0, 0, -g(1) * (t232 * t249 + t243 * t248) - g(2) * (t234 * t249 + t241 * t248) - g(3) * (t239 * t249 + t248 * t278) -g(1) * (-t232 * t248 + t243 * t249) - g(2) * (-t234 * t248 + t241 * t249) - g(3) * (-t239 * t248 + t249 * t278);];
U_reg  = t1;
