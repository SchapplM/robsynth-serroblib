% Calculate minimal parameter regressor of potential energy for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:49:46
% EndTime: 2019-03-09 21:49:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (245->68), mult. (587->94), div. (0->0), fcn. (748->10), ass. (0->43)
t282 = pkin(5) + pkin(10);
t257 = sin(pkin(6));
t261 = sin(qJ(2));
t281 = t257 * t261;
t262 = sin(qJ(1));
t280 = t257 * t262;
t264 = cos(qJ(3));
t279 = t257 * t264;
t265 = cos(qJ(2));
t278 = t257 * t265;
t266 = cos(qJ(1));
t277 = t257 * t266;
t276 = t262 * t261;
t275 = t262 * t265;
t274 = t266 * t261;
t273 = t266 * t265;
t258 = cos(pkin(6));
t245 = t258 * t274 + t275;
t260 = sin(qJ(3));
t234 = t245 * t264 - t260 * t277;
t244 = -t258 * t273 + t276;
t259 = sin(qJ(4));
t263 = cos(qJ(4));
t224 = t234 * t259 - t244 * t263;
t247 = -t258 * t276 + t273;
t236 = t247 * t264 + t260 * t280;
t246 = t258 * t275 + t274;
t226 = t236 * t259 - t246 * t263;
t243 = t258 * t260 + t261 * t279;
t231 = t243 * t259 + t263 * t278;
t272 = g(1) * t226 + g(2) * t224 + g(3) * t231;
t225 = t234 * t263 + t244 * t259;
t227 = t236 * t263 + t246 * t259;
t232 = t243 * t263 - t259 * t278;
t271 = g(1) * t227 + g(2) * t225 + g(3) * t232;
t233 = t245 * t260 + t264 * t277;
t235 = t247 * t260 - t262 * t279;
t242 = -t258 * t264 + t260 * t281;
t270 = g(1) * t235 + g(2) * t233 + g(3) * t242;
t269 = t266 * pkin(1) + t247 * pkin(2) + t236 * pkin(3) + t227 * pkin(4) + pkin(8) * t280 + t246 * pkin(9) + t226 * qJ(5);
t268 = pkin(2) * t281 + t243 * pkin(3) + t232 * pkin(4) + t258 * pkin(8) - pkin(9) * t278 + t231 * qJ(5) + pkin(7);
t267 = t262 * pkin(1) + t245 * pkin(2) + t234 * pkin(3) + t225 * pkin(4) - pkin(8) * t277 + t244 * pkin(9) + t224 * qJ(5);
t1 = [0, -g(1) * t266 - g(2) * t262, g(1) * t262 - g(2) * t266, 0, 0, 0, 0, 0, -g(1) * t247 - g(2) * t245 - g(3) * t281, g(1) * t246 + g(2) * t244 - g(3) * t278, 0, 0, 0, 0, 0, -g(1) * t236 - g(2) * t234 - g(3) * t243, t270, 0, 0, 0, 0, 0, -t271, t272, -t270, t271, -t272, -g(1) * (t235 * pkin(10) + t269) - g(2) * (t233 * pkin(10) + t267) - g(3) * (t242 * pkin(10) + t268) -t270, -t272, -t271, -g(1) * (t227 * qJ(6) + t282 * t235 + t269) - g(2) * (t225 * qJ(6) + t282 * t233 + t267) - g(3) * (t232 * qJ(6) + t282 * t242 + t268);];
U_reg  = t1;
