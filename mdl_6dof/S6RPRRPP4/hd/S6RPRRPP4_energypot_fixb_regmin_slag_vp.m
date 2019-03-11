% Calculate minimal parameter regressor of potential energy for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:40:45
% EndTime: 2019-03-09 04:40:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (153->55), mult. (157->73), div. (0->0), fcn. (160->10), ass. (0->36)
t231 = cos(qJ(4));
t217 = t231 * pkin(4) + pkin(3);
t223 = pkin(9) + qJ(3);
t218 = sin(t223);
t220 = cos(t223);
t226 = cos(pkin(9));
t227 = -qJ(5) - pkin(8);
t248 = t226 * pkin(2) + t217 * t220 - t218 * t227 + pkin(1);
t247 = g(3) * t218;
t224 = qJ(4) + pkin(10);
t219 = sin(t224);
t230 = sin(qJ(1));
t244 = t230 * t219;
t221 = cos(t224);
t243 = t230 * t221;
t229 = sin(qJ(4));
t242 = t230 * t229;
t241 = t230 * t231;
t232 = cos(qJ(1));
t240 = t232 * t219;
t239 = t232 * t221;
t238 = t232 * t229;
t237 = t232 * t231;
t225 = sin(pkin(9));
t236 = t225 * pkin(2) + t218 * t217 + t220 * t227 + pkin(6);
t235 = g(1) * t232 + g(2) * t230;
t228 = -pkin(7) - qJ(2);
t234 = pkin(4) * t242 - t230 * t228 + t248 * t232;
t233 = -pkin(4) * t238 + t232 * t228 + t248 * t230;
t211 = g(1) * t230 - g(2) * t232;
t205 = t220 * t239 + t244;
t204 = t220 * t240 - t243;
t203 = t220 * t243 - t240;
t202 = t220 * t244 + t239;
t201 = -g(3) * t220 + t235 * t218;
t1 = [0, -t235, t211, -g(3) * t225 - t235 * t226, -g(3) * t226 + t235 * t225, -t211, -g(1) * (t232 * pkin(1) + t230 * qJ(2)) - g(2) * (t230 * pkin(1) - t232 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t235 * t220 - t247, t201, 0, 0, 0, 0, 0, -g(1) * (t220 * t237 + t242) - g(2) * (t220 * t241 - t238) - t231 * t247, -g(1) * (-t220 * t238 + t241) - g(2) * (-t220 * t242 - t237) + t229 * t247, -t201, -g(1) * t234 - g(2) * t233 - g(3) * t236, -g(1) * t205 - g(2) * t203 - t221 * t247, -t201, -g(1) * t204 - g(2) * t202 - t219 * t247, -g(1) * (t205 * pkin(5) + t204 * qJ(6) + t234) - g(2) * (t203 * pkin(5) + t202 * qJ(6) + t233) - g(3) * ((pkin(5) * t221 + qJ(6) * t219) * t218 + t236);];
U_reg  = t1;
