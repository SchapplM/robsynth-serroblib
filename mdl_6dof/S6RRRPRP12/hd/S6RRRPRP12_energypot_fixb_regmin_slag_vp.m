% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:59:16
% EndTime: 2019-03-09 17:59:16
% DurationCPUTime: 0.17s
% Computational Cost: add. (198->67), mult. (470->94), div. (0->0), fcn. (590->10), ass. (0->45)
t280 = pkin(4) + pkin(9);
t254 = sin(pkin(6));
t258 = sin(qJ(2));
t279 = t254 * t258;
t259 = sin(qJ(1));
t278 = t254 * t259;
t261 = cos(qJ(3));
t277 = t254 * t261;
t262 = cos(qJ(2));
t276 = t254 * t262;
t263 = cos(qJ(1));
t275 = t254 * t263;
t274 = t258 * t259;
t273 = t258 * t263;
t272 = t259 * t262;
t271 = t262 * t263;
t255 = cos(pkin(6));
t257 = sin(qJ(3));
t239 = -t255 * t261 + t257 * t279;
t240 = t255 * t257 + t258 * t277;
t270 = pkin(2) * t279 + t240 * pkin(3) + t255 * pkin(8) + qJ(4) * t239 + pkin(7);
t244 = -t255 * t274 + t271;
t234 = t244 * t257 - t259 * t277;
t235 = t244 * t261 + t257 * t278;
t269 = t263 * pkin(1) + t244 * pkin(2) + t235 * pkin(3) + pkin(8) * t278 + qJ(4) * t234;
t242 = t255 * t273 + t272;
t232 = t242 * t257 + t261 * t275;
t241 = -t255 * t271 + t274;
t256 = sin(qJ(5));
t260 = cos(qJ(5));
t222 = -t232 * t260 + t241 * t256;
t243 = t255 * t272 + t273;
t224 = -t234 * t260 + t243 * t256;
t230 = t239 * t260 + t256 * t276;
t268 = g(1) * t224 + g(2) * t222 - g(3) * t230;
t267 = g(1) * t234 + g(2) * t232 + g(3) * t239;
t233 = t242 * t261 - t257 * t275;
t266 = g(1) * t235 + g(2) * t233 + g(3) * t240;
t265 = -g(1) * t243 - g(2) * t241 + g(3) * t276;
t264 = t259 * pkin(1) + t242 * pkin(2) + t233 * pkin(3) - pkin(8) * t275 + t232 * qJ(4);
t231 = t239 * t256 - t260 * t276;
t225 = t234 * t256 + t243 * t260;
t223 = t232 * t256 + t241 * t260;
t220 = -g(1) * t225 - g(2) * t223 - g(3) * t231;
t1 = [0, -g(1) * t263 - g(2) * t259, g(1) * t259 - g(2) * t263, 0, 0, 0, 0, 0, -g(1) * t244 - g(2) * t242 - g(3) * t279, -t265, 0, 0, 0, 0, 0, -t266, t267, t265, t266, -t267, -g(1) * (pkin(9) * t243 + t269) - g(2) * (t241 * pkin(9) + t264) - g(3) * (-pkin(9) * t276 + t270) 0, 0, 0, 0, 0, t220, t268, t220, -t266, -t268, -g(1) * (pkin(5) * t225 + pkin(10) * t235 + qJ(6) * t224 + t280 * t243 + t269) - g(2) * (t223 * pkin(5) + t233 * pkin(10) + t222 * qJ(6) + t280 * t241 + t264) - g(3) * (pkin(5) * t231 + pkin(10) * t240 - t230 * qJ(6) - t280 * t276 + t270);];
U_reg  = t1;
