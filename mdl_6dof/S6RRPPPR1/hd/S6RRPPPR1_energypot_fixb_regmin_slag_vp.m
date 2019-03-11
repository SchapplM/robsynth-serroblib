% Calculate minimal parameter regressor of potential energy for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:17
% EndTime: 2019-03-09 08:08:17
% DurationCPUTime: 0.15s
% Computational Cost: add. (142->49), mult. (182->71), div. (0->0), fcn. (199->10), ass. (0->35)
t237 = qJ(2) + pkin(9);
t234 = sin(t237);
t235 = cos(t237);
t262 = pkin(3) * t235 + qJ(4) * t234;
t260 = g(3) * t234;
t242 = sin(qJ(2));
t259 = t242 * pkin(2) + pkin(6);
t238 = sin(pkin(10));
t243 = sin(qJ(1));
t257 = t243 * t238;
t239 = cos(pkin(10));
t256 = t243 * t239;
t246 = cos(qJ(1));
t255 = t246 * t238;
t254 = t246 * t239;
t245 = cos(qJ(2));
t233 = t245 * pkin(2) + pkin(1);
t240 = -pkin(7) - qJ(3);
t253 = t243 * t233 + t246 * t240;
t252 = t246 * t233 - t243 * t240;
t251 = t243 * t262 + t253;
t250 = g(1) * t246 + g(2) * t243;
t249 = t234 * pkin(3) - t235 * qJ(4) + t259;
t248 = t246 * t262 + t252;
t219 = t235 * t257 + t254;
t221 = t235 * t255 - t256;
t247 = g(1) * t221 + g(2) * t219 + t238 * t260;
t244 = cos(qJ(6));
t241 = sin(qJ(6));
t225 = g(1) * t243 - g(2) * t246;
t222 = t235 * t254 + t257;
t220 = t235 * t256 - t255;
t218 = g(3) * t235 - t250 * t234;
t217 = -g(1) * t222 - g(2) * t220 - t239 * t260;
t1 = [0, -t250, t225, 0, 0, 0, 0, 0, -g(3) * t242 - t250 * t245, -g(3) * t245 + t250 * t242, -t225, -g(1) * t252 - g(2) * t253 - g(3) * t259, t217, t247, t218, -g(1) * t248 - g(2) * t251 - g(3) * t249, t217, t218, -t247, -g(1) * (t222 * pkin(4) + t221 * qJ(5) + t248) - g(2) * (t220 * pkin(4) + t219 * qJ(5) + t251) - g(3) * ((pkin(4) * t239 + qJ(5) * t238) * t234 + t249) 0, 0, 0, 0, 0, -g(1) * (t221 * t241 + t222 * t244) - g(2) * (t219 * t241 + t220 * t244) - (t238 * t241 + t239 * t244) * t260, -g(1) * (t221 * t244 - t222 * t241) - g(2) * (t219 * t244 - t220 * t241) - (t238 * t244 - t239 * t241) * t260;];
U_reg  = t1;
