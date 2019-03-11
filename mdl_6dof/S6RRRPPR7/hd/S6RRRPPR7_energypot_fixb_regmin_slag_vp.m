% Calculate minimal parameter regressor of potential energy for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:01:04
% EndTime: 2019-03-09 16:01:05
% DurationCPUTime: 0.14s
% Computational Cost: add. (122->61), mult. (245->91), div. (0->0), fcn. (285->10), ass. (0->31)
t235 = sin(qJ(2));
t252 = g(3) * t235;
t234 = sin(qJ(3));
t251 = t234 * t235;
t236 = sin(qJ(1));
t250 = t235 * t236;
t237 = cos(qJ(3));
t249 = t235 * t237;
t239 = cos(qJ(1));
t248 = t235 * t239;
t238 = cos(qJ(2));
t247 = t236 * t238;
t246 = t239 * t234;
t245 = t239 * t237;
t244 = t235 * pkin(2) + pkin(3) * t249 + qJ(4) * t251 + pkin(6);
t243 = g(1) * t239 + g(2) * t236;
t216 = -t236 * t237 + t238 * t246;
t217 = t236 * t234 + t238 * t245;
t242 = t217 * pkin(3) + t236 * pkin(7) + pkin(8) * t248 + t216 * qJ(4) + (pkin(2) * t238 + pkin(1)) * t239;
t214 = t234 * t247 + t245;
t241 = g(1) * t216 + g(2) * t214 + g(3) * t251;
t215 = t237 * t247 - t246;
t240 = t236 * pkin(1) + pkin(2) * t247 + t215 * pkin(3) - t239 * pkin(7) + pkin(8) * t250 + t214 * qJ(4);
t233 = cos(pkin(10));
t232 = sin(pkin(10));
t231 = pkin(10) + qJ(6);
t226 = cos(t231);
t225 = sin(t231);
t211 = -g(3) * t238 + t243 * t235;
t210 = -g(1) * t217 - g(2) * t215 - g(3) * t249;
t1 = [0, -t243, g(1) * t236 - g(2) * t239, 0, 0, 0, 0, 0, -t243 * t238 - t252, t211, 0, 0, 0, 0, 0, t210, t241, t210, -t211, -t241, -g(1) * t242 - g(2) * t240 - g(3) * (-t238 * pkin(8) + t244) -g(1) * (t216 * t232 + t217 * t233) - g(2) * (t214 * t232 + t215 * t233) - (t232 * t234 + t233 * t237) * t252, -g(1) * (t216 * t233 - t217 * t232) - g(2) * (t214 * t233 - t215 * t232) - (-t232 * t237 + t233 * t234) * t252, t211, -g(1) * (t217 * pkin(4) - qJ(5) * t248 + t242) - g(2) * (t215 * pkin(4) - qJ(5) * t250 + t240) - g(3) * (pkin(4) * t249 + (-pkin(8) + qJ(5)) * t238 + t244) 0, 0, 0, 0, 0, -g(1) * (t216 * t225 + t217 * t226) - g(2) * (t214 * t225 + t215 * t226) - (t225 * t234 + t226 * t237) * t252, -g(1) * (t216 * t226 - t217 * t225) - g(2) * (t214 * t226 - t215 * t225) - (-t225 * t237 + t226 * t234) * t252;];
U_reg  = t1;
