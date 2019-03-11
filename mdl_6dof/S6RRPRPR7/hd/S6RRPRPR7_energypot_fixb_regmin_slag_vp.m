% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:47:42
% EndTime: 2019-03-09 10:47:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (81->45), mult. (142->65), div. (0->0), fcn. (152->10), ass. (0->31)
t234 = sin(qJ(2));
t254 = qJ(3) * t234 + pkin(1);
t235 = sin(qJ(1));
t239 = cos(qJ(1));
t243 = g(1) * t239 + g(2) * t235;
t230 = qJ(4) + pkin(10);
t224 = sin(t230);
t225 = cos(t230);
t238 = cos(qJ(2));
t251 = g(3) * (-t238 * t224 + t234 * t225);
t250 = t234 * pkin(2) + pkin(6);
t233 = sin(qJ(4));
t248 = t234 * t233;
t247 = t235 * t238;
t246 = pkin(4) * t248;
t245 = pkin(2) * t247 + t254 * t235;
t244 = t235 * pkin(7) + (pkin(2) * t238 + t254) * t239;
t242 = t224 * t234 + t225 * t238;
t237 = cos(qJ(4));
t241 = t238 * t233 - t234 * t237;
t240 = t238 * t237 + t248;
t236 = cos(qJ(6));
t232 = sin(qJ(6));
t231 = -qJ(5) - pkin(8);
t223 = t237 * pkin(4) + pkin(3);
t218 = g(1) * t235 - g(2) * t239;
t216 = -g(3) * t234 - t243 * t238;
t215 = -g(3) * t238 + t243 * t234;
t214 = t242 * t239;
t213 = t242 * t235;
t1 = [0, -t243, t218, 0, 0, 0, 0, 0, t216, t215, t216, -t218, -t215, -g(1) * t244 - g(2) * (-t239 * pkin(7) + t245) - g(3) * (-t238 * qJ(3) + t250) 0, 0, 0, 0, 0, g(3) * t241 - t243 * t240, g(3) * t240 + t243 * t241, t218, -g(1) * (t235 * t231 + t244) - g(2) * (t223 * t247 + t235 * t246 + t245) - g(3) * (t234 * t223 + (-pkin(4) * t233 - qJ(3)) * t238 + t250) + (-g(1) * (t223 * t238 + t246) - g(2) * (-pkin(7) - t231)) * t239, 0, 0, 0, 0, 0, -g(1) * (t214 * t236 - t235 * t232) - g(2) * (t213 * t236 + t239 * t232) - t236 * t251, -g(1) * (-t214 * t232 - t235 * t236) - g(2) * (-t213 * t232 + t239 * t236) + t232 * t251;];
U_reg  = t1;
