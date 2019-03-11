% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR14_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:20:35
% EndTime: 2019-03-09 20:20:35
% DurationCPUTime: 0.14s
% Computational Cost: add. (133->61), mult. (296->98), div. (0->0), fcn. (373->12), ass. (0->36)
t256 = sin(pkin(6));
t260 = sin(qJ(2));
t277 = t256 * t260;
t261 = sin(qJ(1));
t276 = t256 * t261;
t263 = cos(qJ(3));
t275 = t256 * t263;
t264 = cos(qJ(2));
t274 = t256 * t264;
t265 = cos(qJ(1));
t273 = t256 * t265;
t272 = t261 * t260;
t271 = t261 * t264;
t270 = t265 * t260;
t269 = t265 * t264;
t257 = cos(pkin(6));
t247 = t257 * t270 + t271;
t259 = sin(qJ(3));
t240 = t247 * t259 + t263 * t273;
t249 = -t257 * t272 + t269;
t242 = t249 * t259 - t261 * t275;
t244 = -t257 * t263 + t259 * t277;
t268 = g(1) * t242 + g(2) * t240 + g(3) * t244;
t241 = t247 * t263 - t259 * t273;
t243 = t249 * t263 + t259 * t276;
t245 = t257 * t259 + t260 * t275;
t267 = g(1) * t243 + g(2) * t241 + g(3) * t245;
t246 = -t257 * t269 + t272;
t248 = t257 * t271 + t270;
t266 = -g(1) * t248 - g(2) * t246 + g(3) * t274;
t262 = cos(qJ(5));
t258 = sin(qJ(5));
t255 = qJ(5) + qJ(6);
t254 = cos(t255);
t253 = sin(t255);
t1 = [0, -g(1) * t265 - g(2) * t261, g(1) * t261 - g(2) * t265, 0, 0, 0, 0, 0, -g(1) * t249 - g(2) * t247 - g(3) * t277, -t266, 0, 0, 0, 0, 0, -t267, t268, t266, t267, -t268, -g(1) * (t265 * pkin(1) + t249 * pkin(2) + t243 * pkin(3) + pkin(8) * t276 + t248 * pkin(9) + t242 * qJ(4)) - g(2) * (t261 * pkin(1) + t247 * pkin(2) + t241 * pkin(3) - pkin(8) * t273 + t246 * pkin(9) + t240 * qJ(4)) - g(3) * (t245 * pkin(3) + t257 * pkin(8) + t244 * qJ(4) + pkin(7) + (pkin(2) * t260 - pkin(9) * t264) * t256) 0, 0, 0, 0, 0, -g(1) * (t242 * t258 + t248 * t262) - g(2) * (t240 * t258 + t246 * t262) - g(3) * (t244 * t258 - t262 * t274) -g(1) * (t242 * t262 - t248 * t258) - g(2) * (t240 * t262 - t246 * t258) - g(3) * (t244 * t262 + t258 * t274) 0, 0, 0, 0, 0, -g(1) * (t242 * t253 + t248 * t254) - g(2) * (t240 * t253 + t246 * t254) - g(3) * (t244 * t253 - t254 * t274) -g(1) * (t242 * t254 - t248 * t253) - g(2) * (t240 * t254 - t246 * t253) - g(3) * (t244 * t254 + t253 * t274);];
U_reg  = t1;
