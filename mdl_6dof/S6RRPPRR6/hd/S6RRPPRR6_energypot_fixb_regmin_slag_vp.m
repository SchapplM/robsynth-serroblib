% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:15
% EndTime: 2019-03-09 09:15:15
% DurationCPUTime: 0.11s
% Computational Cost: add. (95->42), mult. (157->63), div. (0->0), fcn. (173->10), ass. (0->29)
t234 = sin(qJ(2));
t253 = qJ(3) * t234 + pkin(1);
t235 = sin(qJ(1));
t238 = cos(qJ(1));
t243 = g(1) * t238 + g(2) * t235;
t230 = pkin(10) + qJ(5);
t224 = sin(t230);
t225 = cos(t230);
t237 = cos(qJ(2));
t242 = t237 * t224 - t234 * t225;
t250 = g(3) * t242;
t248 = t235 * t237;
t247 = t237 * t238;
t246 = pkin(2) * t248 + t253 * t235;
t245 = pkin(2) * t247 + t235 * pkin(7) + t253 * t238;
t244 = t234 * pkin(2) - t237 * qJ(3) + pkin(6);
t241 = t234 * t224 + t237 * t225;
t231 = sin(pkin(10));
t232 = cos(pkin(10));
t240 = t237 * t231 - t234 * t232;
t239 = t234 * t231 + t237 * t232;
t236 = cos(qJ(6));
t233 = sin(qJ(6));
t219 = g(1) * t235 - g(2) * t238;
t217 = -g(3) * t234 - t243 * t237;
t216 = -g(3) * t237 + t243 * t234;
t215 = t241 * t238;
t214 = t241 * t235;
t1 = [0, -t243, t219, 0, 0, 0, 0, 0, t217, t216, t217, -t219, -t216, -g(1) * t245 - g(2) * (-t238 * pkin(7) + t246) - g(3) * t244, g(3) * t240 - t243 * t239, g(3) * t239 + t243 * t240, t219, -g(1) * (pkin(3) * t247 - t235 * qJ(4) + t245) - g(2) * (pkin(3) * t248 + (-pkin(7) + qJ(4)) * t238 + t246) - g(3) * (t234 * pkin(3) + t244) 0, 0, 0, 0, 0, -g(1) * t215 - g(2) * t214 + t250, g(3) * t241 + t243 * t242, 0, 0, 0, 0, 0, -g(1) * (t215 * t236 - t235 * t233) - g(2) * (t214 * t236 + t238 * t233) + t236 * t250, -g(1) * (-t215 * t233 - t235 * t236) - g(2) * (-t214 * t233 + t238 * t236) - t233 * t250;];
U_reg  = t1;
