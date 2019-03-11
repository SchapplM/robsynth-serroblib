% Calculate minimal parameter regressor of potential energy for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:20:47
% EndTime: 2019-03-08 22:20:47
% DurationCPUTime: 0.14s
% Computational Cost: add. (169->69), mult. (327->115), div. (0->0), fcn. (416->14), ass. (0->35)
t241 = sin(pkin(6));
t256 = pkin(7) * t241;
t245 = sin(qJ(3));
t255 = t241 * t245;
t246 = sin(qJ(2));
t254 = t241 * t246;
t247 = cos(qJ(3));
t253 = t241 * t247;
t248 = cos(qJ(2));
t252 = t241 * t248;
t244 = cos(pkin(6));
t251 = t244 * t246;
t250 = t244 * t248;
t238 = pkin(12) + qJ(5);
t240 = sin(pkin(11));
t243 = cos(pkin(11));
t225 = t240 * t248 + t243 * t251;
t220 = t225 * t245 + t243 * t253;
t227 = -t240 * t251 + t243 * t248;
t222 = t227 * t245 - t240 * t253;
t228 = -t244 * t247 + t245 * t254;
t249 = g(1) * t222 + g(2) * t220 + g(3) * t228;
t242 = cos(pkin(12));
t239 = sin(pkin(12));
t237 = qJ(6) + t238;
t236 = cos(t238);
t235 = sin(t238);
t234 = cos(t237);
t233 = sin(t237);
t229 = t244 * t245 + t246 * t253;
t226 = t240 * t250 + t243 * t246;
t224 = t240 * t246 - t243 * t250;
t223 = t227 * t247 + t240 * t255;
t221 = t225 * t247 - t243 * t255;
t1 = [-g(3) * qJ(1), 0, -g(1) * t227 - g(2) * t225 - g(3) * t254, g(1) * t226 + g(2) * t224 - g(3) * t252, 0, 0, 0, 0, 0, -g(1) * t223 - g(2) * t221 - g(3) * t229, t249, -g(1) * (t223 * t242 + t226 * t239) - g(2) * (t221 * t242 + t224 * t239) - g(3) * (t229 * t242 - t239 * t252) -g(1) * (-t223 * t239 + t226 * t242) - g(2) * (-t221 * t239 + t224 * t242) - g(3) * (-t229 * t239 - t242 * t252) -t249, -g(1) * (t243 * pkin(1) + t227 * pkin(2) + t223 * pkin(3) + t226 * pkin(8) + t222 * qJ(4) + t240 * t256) - g(2) * (t240 * pkin(1) + t225 * pkin(2) + t221 * pkin(3) + t224 * pkin(8) + t220 * qJ(4) - t243 * t256) - g(3) * (t229 * pkin(3) + t244 * pkin(7) + t228 * qJ(4) + qJ(1) + (pkin(2) * t246 - pkin(8) * t248) * t241) 0, 0, 0, 0, 0, -g(1) * (t223 * t236 + t226 * t235) - g(2) * (t221 * t236 + t224 * t235) - g(3) * (t229 * t236 - t235 * t252) -g(1) * (-t223 * t235 + t226 * t236) - g(2) * (-t221 * t235 + t224 * t236) - g(3) * (-t229 * t235 - t236 * t252) 0, 0, 0, 0, 0, -g(1) * (t223 * t234 + t226 * t233) - g(2) * (t221 * t234 + t224 * t233) - g(3) * (t229 * t234 - t233 * t252) -g(1) * (-t223 * t233 + t226 * t234) - g(2) * (-t221 * t233 + t224 * t234) - g(3) * (-t229 * t233 - t234 * t252);];
U_reg  = t1;
