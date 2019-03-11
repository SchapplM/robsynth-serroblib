% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP14_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:09:50
% EndTime: 2019-03-09 13:09:50
% DurationCPUTime: 0.15s
% Computational Cost: add. (169->68), mult. (400->98), div. (0->0), fcn. (497->10), ass. (0->46)
t254 = sin(qJ(1));
t276 = g(1) * t254;
t258 = cos(qJ(1));
t275 = g(2) * t258;
t249 = sin(pkin(6));
t253 = sin(qJ(2));
t274 = t249 * t253;
t273 = t249 * t254;
t257 = cos(qJ(2));
t272 = t249 * t257;
t271 = t249 * t258;
t270 = t254 * t253;
t269 = t254 * t257;
t268 = t258 * t253;
t267 = t258 * t257;
t250 = cos(pkin(6));
t266 = pkin(2) * t274 + t250 * pkin(8) + pkin(7);
t265 = -t275 + t276;
t237 = -t250 * t267 + t270;
t238 = t250 * t268 + t269;
t264 = t254 * pkin(1) + t238 * pkin(2) + t237 * qJ(3);
t239 = t250 * t269 + t268;
t240 = -t250 * t270 + t267;
t263 = t258 * pkin(1) + t240 * pkin(2) + pkin(8) * t273 + t239 * qJ(3);
t252 = sin(qJ(4));
t256 = cos(qJ(4));
t227 = t239 * t252 + t256 * t273;
t251 = sin(qJ(5));
t255 = cos(qJ(5));
t220 = t227 * t251 - t240 * t255;
t229 = t237 * t252 - t256 * t271;
t222 = t229 * t251 - t238 * t255;
t236 = t250 * t256 - t252 * t272;
t224 = t236 * t251 - t255 * t274;
t262 = g(1) * t220 + g(2) * t222 + g(3) * t224;
t226 = -t239 * t256 + t252 * t273;
t228 = t237 * t256 + t252 * t271;
t235 = t250 * t252 + t256 * t272;
t261 = g(1) * t226 - g(2) * t228 + g(3) * t235;
t260 = -g(1) * t239 - g(2) * t237 + g(3) * t272;
t259 = g(1) * t240 + g(2) * t238 + g(3) * t274;
t225 = t236 * t255 + t251 * t274;
t223 = t229 * t255 + t238 * t251;
t221 = t227 * t255 + t240 * t251;
t219 = -g(1) * t221 - g(2) * t223 - g(3) * t225;
t1 = [0, -g(1) * t258 - g(2) * t254, t265, 0, 0, 0, 0, 0, -t259, -t260, -g(3) * t250 - t265 * t249, t259, t260, -g(1) * t263 - g(2) * (-pkin(8) * t271 + t264) - g(3) * (-qJ(3) * t272 + t266) 0, 0, 0, 0, 0, -g(1) * t227 - g(2) * t229 - g(3) * t236, t261, 0, 0, 0, 0, 0, t219, t262, t219, -t261, -t262, -g(1) * (t227 * pkin(4) + t221 * pkin(5) + t240 * pkin(9) + t226 * pkin(10) + t220 * qJ(6) + t263) - g(2) * (t229 * pkin(4) + t223 * pkin(5) + t238 * pkin(9) - t228 * pkin(10) + t222 * qJ(6) + t264) - g(3) * (t250 * pkin(3) + t236 * pkin(4) + t225 * pkin(5) + t235 * pkin(10) + t224 * qJ(6) + t266) + (-pkin(3) * t276 - g(3) * (pkin(9) * t253 - qJ(3) * t257) - (-pkin(3) - pkin(8)) * t275) * t249;];
U_reg  = t1;
