% Calculate minimal parameter regressor of potential energy for
% S6PRRPRR2
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
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:16
% EndTime: 2019-03-08 22:00:16
% DurationCPUTime: 0.14s
% Computational Cost: add. (135->60), mult. (232->102), div. (0->0), fcn. (288->14), ass. (0->35)
t234 = qJ(3) + pkin(12);
t230 = sin(t234);
t237 = sin(pkin(6));
t256 = t230 * t237;
t242 = sin(qJ(3));
t255 = t237 * t242;
t243 = sin(qJ(2));
t254 = t237 * t243;
t245 = cos(qJ(3));
t253 = t237 * t245;
t246 = cos(qJ(2));
t252 = t237 * t246;
t239 = cos(pkin(6));
t251 = t239 * t242;
t250 = t239 * t243;
t249 = t239 * t246;
t236 = sin(pkin(11));
t238 = cos(pkin(11));
t224 = t236 * t243 - t238 * t249;
t226 = t236 * t249 + t238 * t243;
t247 = -g(1) * t226 - g(2) * t224 + g(3) * t252;
t244 = cos(qJ(5));
t241 = sin(qJ(5));
t240 = -qJ(4) - pkin(8);
t235 = qJ(5) + qJ(6);
t233 = cos(t235);
t232 = sin(t235);
t231 = cos(t234);
t229 = t245 * pkin(3) + pkin(2);
t227 = -t236 * t250 + t238 * t246;
t225 = t236 * t246 + t238 * t250;
t223 = t239 * t230 + t231 * t254;
t222 = t227 * t231 + t236 * t256;
t221 = t225 * t231 - t238 * t256;
t1 = [-g(3) * qJ(1), 0, -g(1) * t227 - g(2) * t225 - g(3) * t254, -t247, 0, 0, 0, 0, 0, -g(1) * (t227 * t245 + t236 * t255) - g(2) * (t225 * t245 - t238 * t255) - g(3) * (t243 * t253 + t251) -g(1) * (-t227 * t242 + t236 * t253) - g(2) * (-t225 * t242 - t238 * t253) - g(3) * (t239 * t245 - t242 * t254) t247, -g(1) * (t238 * pkin(1) - t226 * t240 + t227 * t229) - g(2) * (t236 * pkin(1) - t224 * t240 + t225 * t229) - g(3) * (pkin(3) * t251 + t239 * pkin(7) + qJ(1)) + (-g(3) * (t229 * t243 + t240 * t246) + (-g(1) * t236 + g(2) * t238) * (pkin(3) * t242 + pkin(7))) * t237, 0, 0, 0, 0, 0, -g(1) * (t222 * t244 + t226 * t241) - g(2) * (t221 * t244 + t224 * t241) - g(3) * (t223 * t244 - t241 * t252) -g(1) * (-t222 * t241 + t226 * t244) - g(2) * (-t221 * t241 + t224 * t244) - g(3) * (-t223 * t241 - t244 * t252) 0, 0, 0, 0, 0, -g(1) * (t222 * t233 + t226 * t232) - g(2) * (t221 * t233 + t224 * t232) - g(3) * (t223 * t233 - t232 * t252) -g(1) * (-t222 * t232 + t226 * t233) - g(2) * (-t221 * t232 + t224 * t233) - g(3) * (-t223 * t232 - t233 * t252);];
U_reg  = t1;
