% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:11:38
% EndTime: 2019-03-10 02:11:39
% DurationCPUTime: 0.18s
% Computational Cost: add. (148->63), mult. (290->102), div. (0->0), fcn. (362->12), ass. (0->37)
t250 = qJ(4) + qJ(5);
t247 = sin(t250);
t253 = sin(qJ(4));
t271 = t253 * pkin(4) + pkin(5) * t247 + pkin(9);
t251 = sin(pkin(6));
t255 = sin(qJ(2));
t270 = t251 * t255;
t258 = cos(qJ(3));
t269 = t251 * t258;
t259 = cos(qJ(2));
t268 = t251 * t259;
t260 = cos(qJ(1));
t267 = t251 * t260;
t256 = sin(qJ(1));
t266 = t256 * t255;
t265 = t256 * t259;
t264 = t260 * t255;
t263 = t260 * t259;
t262 = g(1) * t256 - g(2) * t260;
t252 = cos(pkin(6));
t239 = t252 * t264 + t265;
t254 = sin(qJ(3));
t232 = t239 * t254 + t258 * t267;
t241 = -t252 * t266 + t263;
t234 = t241 * t254 - t256 * t269;
t236 = -t252 * t258 + t254 * t270;
t261 = g(1) * t234 + g(2) * t232 + g(3) * t236;
t257 = cos(qJ(4));
t249 = -qJ(6) - pkin(11) - pkin(10);
t248 = cos(t250);
t242 = t257 * pkin(4) + pkin(5) * t248 + pkin(3);
t240 = t252 * t265 + t264;
t238 = -t252 * t263 + t266;
t237 = t252 * t254 + t255 * t269;
t235 = t256 * t251 * t254 + t241 * t258;
t233 = t239 * t258 - t254 * t267;
t1 = [0, -g(1) * t260 - g(2) * t256, t262, 0, 0, 0, 0, 0, -g(1) * t241 - g(2) * t239 - g(3) * t270, g(1) * t240 + g(2) * t238 - g(3) * t268, 0, 0, 0, 0, 0, -g(1) * t235 - g(2) * t233 - g(3) * t237, t261, 0, 0, 0, 0, 0, -g(1) * (t235 * t257 + t240 * t253) - g(2) * (t233 * t257 + t238 * t253) - g(3) * (t237 * t257 - t253 * t268) -g(1) * (-t235 * t253 + t240 * t257) - g(2) * (-t233 * t253 + t238 * t257) - g(3) * (-t237 * t253 - t257 * t268) 0, 0, 0, 0, 0, -g(1) * (t235 * t248 + t240 * t247) - g(2) * (t233 * t248 + t238 * t247) - g(3) * (t237 * t248 - t247 * t268) -g(1) * (-t235 * t247 + t240 * t248) - g(2) * (-t233 * t247 + t238 * t248) - g(3) * (-t237 * t247 - t248 * t268) -t261, -g(1) * (t260 * pkin(1) + t241 * pkin(2) - t234 * t249 + t235 * t242 + t271 * t240) - g(2) * (t256 * pkin(1) + t239 * pkin(2) - t232 * t249 + t233 * t242 + t271 * t238) - g(3) * (t252 * pkin(8) - t236 * t249 + t237 * t242 + pkin(7)) + (-g(3) * (pkin(2) * t255 - t271 * t259) - t262 * pkin(8)) * t251;];
U_reg  = t1;
