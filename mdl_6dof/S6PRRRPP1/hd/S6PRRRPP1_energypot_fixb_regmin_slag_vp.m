% Calculate minimal parameter regressor of potential energy for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:46:53
% EndTime: 2019-03-08 22:46:53
% DurationCPUTime: 0.18s
% Computational Cost: add. (219->72), mult. (443->109), div. (0->0), fcn. (552->12), ass. (0->45)
t253 = sin(pkin(6));
t275 = pkin(7) * t253;
t252 = sin(pkin(10));
t254 = cos(pkin(10));
t259 = sin(qJ(2));
t255 = cos(pkin(6));
t262 = cos(qJ(2));
t267 = t255 * t262;
t234 = t252 * t259 - t254 * t267;
t257 = sin(qJ(4));
t274 = t234 * t257;
t236 = t252 * t267 + t254 * t259;
t273 = t236 * t257;
t258 = sin(qJ(3));
t272 = t253 * t258;
t271 = t253 * t259;
t261 = cos(qJ(3));
t270 = t253 * t261;
t269 = t253 * t262;
t268 = t255 * t259;
t235 = t252 * t262 + t254 * t268;
t224 = t235 * t258 + t254 * t270;
t237 = -t252 * t268 + t254 * t262;
t226 = t237 * t258 - t252 * t270;
t238 = -t255 * t261 + t258 * t271;
t266 = g(1) * t226 + g(2) * t224 + g(3) * t238;
t227 = t237 * t261 + t252 * t272;
t260 = cos(qJ(4));
t245 = t260 * pkin(4) + pkin(3);
t256 = -qJ(5) - pkin(9);
t265 = t254 * pkin(1) + t237 * pkin(2) + pkin(4) * t273 + t236 * pkin(8) - t226 * t256 + t227 * t245 + t252 * t275;
t225 = t235 * t261 - t254 * t272;
t264 = t252 * pkin(1) + t235 * pkin(2) + pkin(4) * t274 + t234 * pkin(8) - t224 * t256 + t225 * t245 - t254 * t275;
t239 = t255 * t258 + t259 * t270;
t263 = qJ(1) + t239 * t245 - t238 * t256 + pkin(2) * t271 + t255 * pkin(7) + (-pkin(4) * t257 - pkin(8)) * t269;
t251 = qJ(4) + pkin(11);
t247 = cos(t251);
t246 = sin(t251);
t221 = t239 * t247 - t246 * t269;
t220 = t239 * t246 + t247 * t269;
t217 = t227 * t247 + t236 * t246;
t216 = t227 * t246 - t236 * t247;
t215 = t225 * t247 + t234 * t246;
t214 = t225 * t246 - t234 * t247;
t1 = [-g(3) * qJ(1), 0, -g(1) * t237 - g(2) * t235 - g(3) * t271, g(1) * t236 + g(2) * t234 - g(3) * t269, 0, 0, 0, 0, 0, -g(1) * t227 - g(2) * t225 - g(3) * t239, t266, 0, 0, 0, 0, 0, -g(1) * (t227 * t260 + t273) - g(2) * (t225 * t260 + t274) - g(3) * (t239 * t260 - t257 * t269) -g(1) * (-t227 * t257 + t236 * t260) - g(2) * (-t225 * t257 + t234 * t260) - g(3) * (-t239 * t257 - t260 * t269) -t266, -g(1) * t265 - g(2) * t264 - g(3) * t263, -g(1) * t217 - g(2) * t215 - g(3) * t221, -t266, -g(1) * t216 - g(2) * t214 - g(3) * t220, -g(1) * (pkin(5) * t217 + qJ(6) * t216 + t265) - g(2) * (pkin(5) * t215 + qJ(6) * t214 + t264) - g(3) * (t221 * pkin(5) + t220 * qJ(6) + t263);];
U_reg  = t1;
