% Calculate minimal parameter regressor of potential energy for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:51:32
% EndTime: 2019-03-09 00:51:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (131->51), mult. (243->93), div. (0->0), fcn. (316->14), ass. (0->30)
t237 = sin(pkin(6));
t241 = sin(qJ(3));
t251 = t237 * t241;
t242 = sin(qJ(2));
t250 = t237 * t242;
t244 = cos(qJ(3));
t249 = t237 * t244;
t245 = cos(qJ(2));
t248 = t237 * t245;
t239 = cos(pkin(6));
t247 = t239 * t242;
t246 = t239 * t245;
t235 = qJ(4) + qJ(5);
t243 = cos(qJ(4));
t240 = sin(qJ(4));
t238 = cos(pkin(12));
t236 = sin(pkin(12));
t234 = qJ(6) + t235;
t233 = cos(t235);
t232 = sin(t235);
t231 = cos(t234);
t230 = sin(t234);
t229 = t239 * t241 + t242 * t249;
t228 = -t236 * t247 + t238 * t245;
t227 = t236 * t246 + t238 * t242;
t226 = t236 * t245 + t238 * t247;
t225 = t236 * t242 - t238 * t246;
t224 = t228 * t244 + t236 * t251;
t223 = t226 * t244 - t238 * t251;
t1 = [-g(3) * qJ(1), 0, -g(1) * t228 - g(2) * t226 - g(3) * t250, g(1) * t227 + g(2) * t225 - g(3) * t248, 0, 0, 0, 0, 0, -g(1) * t224 - g(2) * t223 - g(3) * t229, -g(1) * (-t228 * t241 + t236 * t249) - g(2) * (-t226 * t241 - t238 * t249) - g(3) * (t239 * t244 - t241 * t250) 0, 0, 0, 0, 0, -g(1) * (t224 * t243 + t227 * t240) - g(2) * (t223 * t243 + t225 * t240) - g(3) * (t229 * t243 - t240 * t248) -g(1) * (-t224 * t240 + t227 * t243) - g(2) * (-t223 * t240 + t225 * t243) - g(3) * (-t229 * t240 - t243 * t248) 0, 0, 0, 0, 0, -g(1) * (t224 * t233 + t227 * t232) - g(2) * (t223 * t233 + t225 * t232) - g(3) * (t229 * t233 - t232 * t248) -g(1) * (-t224 * t232 + t227 * t233) - g(2) * (-t223 * t232 + t225 * t233) - g(3) * (-t229 * t232 - t233 * t248) 0, 0, 0, 0, 0, -g(1) * (t224 * t231 + t227 * t230) - g(2) * (t223 * t231 + t225 * t230) - g(3) * (t229 * t231 - t230 * t248) -g(1) * (-t224 * t230 + t227 * t231) - g(2) * (-t223 * t230 + t225 * t231) - g(3) * (-t229 * t230 - t231 * t248);];
U_reg  = t1;
