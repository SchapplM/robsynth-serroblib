% Calculate minimal parameter regressor of potential energy for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:17:46
% EndTime: 2019-03-09 00:17:46
% DurationCPUTime: 0.12s
% Computational Cost: add. (206->70), mult. (405->108), div. (0->0), fcn. (516->12), ass. (0->43)
t237 = sin(pkin(11));
t239 = cos(pkin(11));
t243 = sin(qJ(2));
t240 = cos(pkin(6));
t246 = cos(qJ(2));
t250 = t240 * t246;
t224 = t237 * t243 - t239 * t250;
t241 = sin(qJ(4));
t257 = t224 * t241;
t226 = t237 * t250 + t239 * t243;
t256 = t226 * t241;
t238 = sin(pkin(6));
t242 = sin(qJ(3));
t255 = t238 * t242;
t254 = t238 * t243;
t245 = cos(qJ(3));
t253 = t238 * t245;
t252 = t238 * t246;
t251 = t240 * t243;
t225 = t237 * t246 + t239 * t251;
t219 = t225 * t245 - t239 * t255;
t236 = qJ(4) + qJ(5);
t234 = sin(t236);
t235 = cos(t236);
t212 = t219 * t234 - t224 * t235;
t227 = -t237 * t251 + t239 * t246;
t221 = t227 * t245 + t237 * t255;
t214 = t221 * t234 - t226 * t235;
t229 = t240 * t242 + t243 * t253;
t216 = t229 * t234 + t235 * t252;
t249 = g(1) * t214 + g(2) * t212 + g(3) * t216;
t218 = t225 * t242 + t239 * t253;
t220 = t227 * t242 - t237 * t253;
t228 = -t240 * t245 + t242 * t254;
t248 = g(1) * t220 + g(2) * t218 + g(3) * t228;
t247 = -pkin(10) - pkin(9);
t244 = cos(qJ(4));
t233 = t244 * pkin(4) + pkin(3);
t217 = t229 * t235 - t234 * t252;
t215 = t221 * t235 + t226 * t234;
t213 = t219 * t235 + t224 * t234;
t211 = -g(1) * t215 - g(2) * t213 - g(3) * t217;
t1 = [-g(3) * qJ(1), 0, -g(1) * t227 - g(2) * t225 - g(3) * t254, g(1) * t226 + g(2) * t224 - g(3) * t252, 0, 0, 0, 0, 0, -g(1) * t221 - g(2) * t219 - g(3) * t229, t248, 0, 0, 0, 0, 0, -g(1) * (t221 * t244 + t256) - g(2) * (t219 * t244 + t257) - g(3) * (t229 * t244 - t241 * t252) -g(1) * (-t221 * t241 + t226 * t244) - g(2) * (-t219 * t241 + t224 * t244) - g(3) * (-t229 * t241 - t244 * t252) 0, 0, 0, 0, 0, t211, t249, t211, -t248, -t249, -g(1) * (t239 * pkin(1) + t227 * pkin(2) + pkin(4) * t256 + t215 * pkin(5) + t226 * pkin(8) + t214 * qJ(6) - t220 * t247 + t221 * t233) - g(2) * (t237 * pkin(1) + t225 * pkin(2) + pkin(4) * t257 + t213 * pkin(5) + t224 * pkin(8) + t212 * qJ(6) - t218 * t247 + t219 * t233) - g(3) * (t217 * pkin(5) + t240 * pkin(7) + t216 * qJ(6) - t228 * t247 + t229 * t233 + qJ(1)) + (-g(3) * (pkin(2) * t243 + (-pkin(4) * t241 - pkin(8)) * t246) + (-g(1) * t237 + g(2) * t239) * pkin(7)) * t238;];
U_reg  = t1;
