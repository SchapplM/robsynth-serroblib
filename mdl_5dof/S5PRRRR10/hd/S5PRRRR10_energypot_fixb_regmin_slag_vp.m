% Calculate minimal parameter regressor of potential energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:05
% EndTime: 2019-12-05 17:26:05
% DurationCPUTime: 0.13s
% Computational Cost: add. (139->46), mult. (383->91), div. (0->0), fcn. (504->14), ass. (0->40)
t230 = sin(pkin(6));
t231 = sin(pkin(5));
t251 = t230 * t231;
t234 = cos(pkin(5));
t250 = t230 * t234;
t233 = cos(pkin(6));
t249 = t231 * t233;
t242 = cos(qJ(2));
t248 = t231 * t242;
t247 = t233 * t242;
t238 = sin(qJ(2));
t246 = t234 * t238;
t245 = t234 * t242;
t229 = sin(pkin(11));
t232 = cos(pkin(11));
t225 = -t229 * t238 + t232 * t245;
t244 = -t225 * t233 + t232 * t251;
t227 = -t229 * t245 - t232 * t238;
t243 = t227 * t233 + t229 * t251;
t241 = cos(qJ(3));
t240 = cos(qJ(4));
t239 = cos(qJ(5));
t237 = sin(qJ(3));
t236 = sin(qJ(4));
t235 = sin(qJ(5));
t228 = -t229 * t246 + t232 * t242;
t226 = t229 * t242 + t232 * t246;
t224 = -t230 * t248 + t234 * t233;
t223 = -t227 * t230 + t229 * t249;
t222 = -t225 * t230 - t232 * t249;
t221 = t237 * t250 + (t237 * t247 + t238 * t241) * t231;
t220 = -t241 * t250 + (t237 * t238 - t241 * t247) * t231;
t219 = t221 * t240 + t224 * t236;
t218 = t228 * t241 + t237 * t243;
t217 = t228 * t237 - t241 * t243;
t216 = t226 * t241 - t237 * t244;
t215 = t226 * t237 + t241 * t244;
t214 = t218 * t240 + t223 * t236;
t213 = t216 * t240 + t222 * t236;
t1 = [-g(3) * qJ(1), 0, -g(3) * t231 * t238 - g(1) * t228 - g(2) * t226, -g(1) * t227 - g(2) * t225 - g(3) * t248, 0, 0, 0, 0, 0, -g(1) * t218 - g(2) * t216 - g(3) * t221, g(1) * t217 + g(2) * t215 + g(3) * t220, 0, 0, 0, 0, 0, -g(1) * t214 - g(2) * t213 - g(3) * t219, -g(1) * (-t218 * t236 + t223 * t240) - g(2) * (-t216 * t236 + t222 * t240) - g(3) * (-t221 * t236 + t224 * t240), 0, 0, 0, 0, 0, -g(1) * (t214 * t239 + t217 * t235) - g(2) * (t213 * t239 + t215 * t235) - g(3) * (t219 * t239 + t220 * t235), -g(1) * (-t214 * t235 + t217 * t239) - g(2) * (-t213 * t235 + t215 * t239) - g(3) * (-t219 * t235 + t220 * t239);];
U_reg = t1;
