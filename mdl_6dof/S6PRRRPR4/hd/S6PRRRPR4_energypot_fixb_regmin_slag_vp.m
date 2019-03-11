% Calculate minimal parameter regressor of potential energy for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:20:27
% EndTime: 2019-03-08 23:20:28
% DurationCPUTime: 0.13s
% Computational Cost: add. (144->63), mult. (281->102), div. (0->0), fcn. (352->12), ass. (0->35)
t236 = sin(pkin(11));
t238 = cos(pkin(11));
t243 = sin(qJ(2));
t239 = cos(pkin(6));
t246 = cos(qJ(2));
t248 = t239 * t246;
t223 = t236 * t243 - t238 * t248;
t241 = sin(qJ(4));
t255 = t223 * t241;
t225 = t236 * t248 + t238 * t243;
t254 = t225 * t241;
t237 = sin(pkin(6));
t242 = sin(qJ(3));
t253 = t237 * t242;
t252 = t237 * t243;
t245 = cos(qJ(3));
t251 = t237 * t245;
t250 = t237 * t246;
t249 = t239 * t243;
t224 = t236 * t246 + t238 * t249;
t219 = t224 * t242 + t238 * t251;
t226 = -t236 * t249 + t238 * t246;
t221 = t226 * t242 - t236 * t251;
t227 = -t239 * t245 + t242 * t252;
t247 = g(1) * t221 + g(2) * t219 + g(3) * t227;
t244 = cos(qJ(4));
t240 = -qJ(5) - pkin(9);
t235 = qJ(4) + pkin(12) + qJ(6);
t234 = t244 * pkin(4) + pkin(3);
t233 = cos(t235);
t232 = sin(t235);
t228 = t239 * t242 + t243 * t251;
t222 = t226 * t245 + t236 * t253;
t220 = t224 * t245 - t238 * t253;
t1 = [-g(3) * qJ(1), 0, -g(1) * t226 - g(2) * t224 - g(3) * t252, g(1) * t225 + g(2) * t223 - g(3) * t250, 0, 0, 0, 0, 0, -g(1) * t222 - g(2) * t220 - g(3) * t228, t247, 0, 0, 0, 0, 0, -g(1) * (t222 * t244 + t254) - g(2) * (t220 * t244 + t255) - g(3) * (t228 * t244 - t241 * t250) -g(1) * (-t222 * t241 + t225 * t244) - g(2) * (-t220 * t241 + t223 * t244) - g(3) * (-t228 * t241 - t244 * t250) -t247, -g(1) * (t238 * pkin(1) + t226 * pkin(2) + pkin(4) * t254 + t225 * pkin(8) - t221 * t240 + t222 * t234) - g(2) * (t236 * pkin(1) + t224 * pkin(2) + pkin(4) * t255 + t223 * pkin(8) - t219 * t240 + t220 * t234) - g(3) * (t239 * pkin(7) - t227 * t240 + t228 * t234 + qJ(1)) + (-g(3) * (pkin(2) * t243 + (-pkin(4) * t241 - pkin(8)) * t246) + (-g(1) * t236 + g(2) * t238) * pkin(7)) * t237, 0, 0, 0, 0, 0, -g(1) * (t222 * t233 + t225 * t232) - g(2) * (t220 * t233 + t223 * t232) - g(3) * (t228 * t233 - t232 * t250) -g(1) * (-t222 * t232 + t225 * t233) - g(2) * (-t220 * t232 + t223 * t233) - g(3) * (-t228 * t232 - t233 * t250);];
U_reg  = t1;
