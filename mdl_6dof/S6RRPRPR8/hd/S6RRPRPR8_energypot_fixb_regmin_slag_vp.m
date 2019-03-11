% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:53:19
% EndTime: 2019-03-09 10:53:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (140->61), mult. (191->87), div. (0->0), fcn. (215->10), ass. (0->32)
t231 = sin(qJ(2));
t247 = g(3) * t231;
t227 = sin(pkin(10));
t232 = sin(qJ(1));
t246 = t232 * t227;
t234 = cos(qJ(2));
t245 = t232 * t234;
t226 = pkin(10) + qJ(4);
t220 = sin(t226);
t235 = cos(qJ(1));
t244 = t235 * t220;
t221 = cos(t226);
t243 = t235 * t221;
t242 = t235 * t227;
t228 = cos(pkin(10));
t241 = t235 * t228;
t240 = t235 * pkin(1) + t232 * pkin(7);
t239 = t232 * pkin(1) - t235 * pkin(7);
t238 = g(1) * t235 + g(2) * t232;
t237 = pkin(2) * t234 + qJ(3) * t231;
t213 = t220 * t245 + t243;
t215 = -t232 * t221 + t234 * t244;
t236 = g(1) * t215 + g(2) * t213 + t220 * t247;
t233 = cos(qJ(6));
t230 = sin(qJ(6));
t229 = -pkin(8) - qJ(3);
t219 = t228 * pkin(3) + pkin(2);
t217 = -g(3) * t234 + t238 * t231;
t216 = t232 * t220 + t234 * t243;
t214 = t221 * t245 - t244;
t212 = -g(1) * t216 - g(2) * t214 - t221 * t247;
t1 = [0, -t238, g(1) * t232 - g(2) * t235, 0, 0, 0, 0, 0, -t238 * t234 - t247, t217, -g(1) * (t234 * t241 + t246) - g(2) * (t228 * t245 - t242) - t228 * t247, -g(1) * (t232 * t228 - t234 * t242) - g(2) * (-t227 * t245 - t241) + t227 * t247, -t217, -g(1) * (t237 * t235 + t240) - g(2) * (t237 * t232 + t239) - g(3) * (t231 * pkin(2) - t234 * qJ(3) + pkin(6)) 0, 0, 0, 0, 0, t212, t236, t212, -t217, -t236, -g(1) * (t235 * t234 * t219 + pkin(3) * t246 + t216 * pkin(4) + t215 * qJ(5) + t240) - g(2) * (-pkin(3) * t242 + t214 * pkin(4) + t213 * qJ(5) + t219 * t245 + t239) - g(3) * (t234 * t229 + pkin(6)) + (-g(3) * (pkin(4) * t221 + qJ(5) * t220 + t219) + t238 * t229) * t231, 0, 0, 0, 0, 0, -g(1) * (t215 * t230 + t216 * t233) - g(2) * (t213 * t230 + t214 * t233) - (t220 * t230 + t221 * t233) * t247, -g(1) * (t215 * t233 - t216 * t230) - g(2) * (t213 * t233 - t214 * t230) - (t220 * t233 - t221 * t230) * t247;];
U_reg  = t1;
