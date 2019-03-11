% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR11
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:03
% EndTime: 2019-03-09 11:15:03
% DurationCPUTime: 0.15s
% Computational Cost: add. (82->49), mult. (122->64), div. (0->0), fcn. (124->8), ass. (0->32)
t227 = sin(qJ(2));
t249 = qJ(3) * t227 + pkin(1);
t226 = sin(qJ(4));
t248 = pkin(4) * t226;
t230 = cos(qJ(2));
t247 = g(3) * t230;
t246 = t227 * pkin(2) + pkin(6);
t220 = qJ(4) + pkin(10) + qJ(6);
t217 = sin(t220);
t228 = sin(qJ(1));
t244 = t228 * t217;
t218 = cos(t220);
t243 = t228 * t218;
t242 = t228 * t226;
t229 = cos(qJ(4));
t241 = t228 * t229;
t240 = t228 * t230;
t231 = cos(qJ(1));
t239 = t231 * t217;
t238 = t231 * t218;
t237 = t231 * t226;
t236 = t231 * t229;
t235 = t227 * t242;
t234 = pkin(2) * t240 + t249 * t228;
t233 = t228 * pkin(7) + (pkin(2) * t230 + t249) * t231;
t232 = g(1) * t231 + g(2) * t228;
t225 = -qJ(5) - pkin(8);
t219 = t229 * pkin(4) + pkin(3);
t212 = g(1) * t228 - g(2) * t231;
t211 = g(3) * t227 + t232 * t230;
t210 = t232 * t227 - t247;
t1 = [0, -t232, t212, 0, 0, 0, 0, 0, -t211, t210, -t212, t211, -t210, -g(1) * t233 - g(2) * (-t231 * pkin(7) + t234) - g(3) * (-t230 * qJ(3) + t246) 0, 0, 0, 0, 0, -g(1) * (t227 * t237 + t241) - g(2) * (t235 - t236) + t226 * t247, -g(1) * (t227 * t236 - t242) - g(2) * (t227 * t241 + t237) + t229 * t247, -t211, -g(1) * (t228 * t219 + t233) - g(2) * (pkin(4) * t235 - t225 * t240 + t234) - g(3) * (-t227 * t225 + (-qJ(3) - t248) * t230 + t246) + (-g(1) * (-t225 * t230 + t227 * t248) - g(2) * (-pkin(7) - t219)) * t231, 0, 0, 0, 0, 0, -g(1) * (t227 * t239 + t243) - g(2) * (t227 * t244 - t238) + t217 * t247, -g(1) * (t227 * t238 - t244) - g(2) * (t227 * t243 + t239) + t218 * t247;];
U_reg  = t1;
