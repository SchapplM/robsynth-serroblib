% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:32:06
% EndTime: 2019-03-09 09:32:06
% DurationCPUTime: 0.16s
% Computational Cost: add. (106->55), mult. (247->87), div. (0->0), fcn. (293->10), ass. (0->36)
t230 = sin(qJ(1));
t250 = g(1) * t230;
t234 = cos(qJ(1));
t249 = g(2) * t234;
t225 = sin(pkin(6));
t229 = sin(qJ(2));
t248 = t225 * t229;
t247 = t225 * t230;
t232 = cos(qJ(5));
t246 = t225 * t232;
t233 = cos(qJ(2));
t245 = t225 * t233;
t244 = t225 * t234;
t243 = t230 * t229;
t242 = t230 * t233;
t241 = t234 * t229;
t240 = t234 * t233;
t226 = cos(pkin(6));
t239 = pkin(2) * t248 + t226 * pkin(8) + pkin(7);
t238 = -t249 + t250;
t214 = -t226 * t240 + t243;
t215 = t226 * t241 + t242;
t237 = t230 * pkin(1) + t215 * pkin(2) + t214 * qJ(3);
t216 = t226 * t242 + t241;
t217 = -t226 * t243 + t240;
t236 = t234 * pkin(1) + t217 * pkin(2) + pkin(8) * t247 + t216 * qJ(3);
t206 = -g(1) * t216 - g(2) * t214 + g(3) * t245;
t235 = g(1) * t217 + g(2) * t215 + g(3) * t248;
t231 = cos(qJ(6));
t228 = sin(qJ(5));
t227 = sin(qJ(6));
t213 = t226 * t232 + t228 * t248;
t210 = -g(3) * t226 - t238 * t225;
t209 = t215 * t228 - t232 * t244;
t208 = t217 * t228 + t230 * t246;
t1 = [0, -g(1) * t234 - g(2) * t230, t238, 0, 0, 0, 0, 0, -t235, -t206, t210, t235, t206, -g(1) * t236 - g(2) * (-pkin(8) * t244 + t237) - g(3) * (-qJ(3) * t245 + t239) t210, t206, -t235, -g(1) * (t217 * qJ(4) + t236) - g(2) * (t215 * qJ(4) + t237) - g(3) * (t226 * pkin(3) + t239) + (-pkin(3) * t250 - g(3) * (-qJ(3) * t233 + qJ(4) * t229) - (-pkin(3) - pkin(8)) * t249) * t225, 0, 0, 0, 0, 0, -g(1) * t208 - g(2) * t209 - g(3) * t213, -g(1) * (t217 * t232 - t228 * t247) - g(2) * (t215 * t232 + t228 * t244) - g(3) * (-t226 * t228 + t229 * t246) 0, 0, 0, 0, 0, -g(1) * (t208 * t231 - t216 * t227) - g(2) * (t209 * t231 - t214 * t227) - g(3) * (t213 * t231 + t227 * t245) -g(1) * (-t208 * t227 - t216 * t231) - g(2) * (-t209 * t227 - t214 * t231) - g(3) * (-t213 * t227 + t231 * t245);];
U_reg  = t1;
