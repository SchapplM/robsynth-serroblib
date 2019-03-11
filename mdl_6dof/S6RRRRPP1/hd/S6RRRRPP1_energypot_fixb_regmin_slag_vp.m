% Calculate minimal parameter regressor of potential energy for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:46:47
% EndTime: 2019-03-09 20:46:47
% DurationCPUTime: 0.13s
% Computational Cost: add. (148->50), mult. (148->66), div. (0->0), fcn. (154->10), ass. (0->35)
t230 = cos(qJ(4));
t216 = t230 * pkin(4) + pkin(3);
t225 = qJ(2) + qJ(3);
t221 = sin(t225);
t222 = cos(t225);
t226 = -qJ(5) - pkin(9);
t231 = cos(qJ(2));
t249 = t231 * pkin(2) + t216 * t222 - t221 * t226 + pkin(1);
t248 = g(3) * t221;
t224 = qJ(4) + pkin(10);
t219 = sin(t224);
t229 = sin(qJ(1));
t245 = t229 * t219;
t220 = cos(t224);
t244 = t229 * t220;
t227 = sin(qJ(4));
t243 = t229 * t227;
t242 = t229 * t230;
t232 = cos(qJ(1));
t241 = t232 * t219;
t240 = t232 * t220;
t239 = t232 * t227;
t238 = t232 * t230;
t228 = sin(qJ(2));
t237 = t228 * pkin(2) + t221 * t216 + t222 * t226 + pkin(6);
t236 = g(1) * t232 + g(2) * t229;
t233 = -pkin(8) - pkin(7);
t235 = pkin(4) * t243 - t229 * t233 + t249 * t232;
t234 = -pkin(4) * t239 + t249 * t229 + t232 * t233;
t207 = t222 * t240 + t245;
t206 = t222 * t241 - t244;
t205 = t222 * t244 - t241;
t204 = t222 * t245 + t240;
t203 = -g(3) * t222 + t236 * t221;
t1 = [0, -t236, g(1) * t229 - g(2) * t232, 0, 0, 0, 0, 0, -g(3) * t228 - t236 * t231, -g(3) * t231 + t236 * t228, 0, 0, 0, 0, 0, -t236 * t222 - t248, t203, 0, 0, 0, 0, 0, -g(1) * (t222 * t238 + t243) - g(2) * (t222 * t242 - t239) - t230 * t248, -g(1) * (-t222 * t239 + t242) - g(2) * (-t222 * t243 - t238) + t227 * t248, -t203, -g(1) * t235 - g(2) * t234 - g(3) * t237, -g(1) * t207 - g(2) * t205 - t220 * t248, -t203, -g(1) * t206 - g(2) * t204 - t219 * t248, -g(1) * (t207 * pkin(5) + t206 * qJ(6) + t235) - g(2) * (t205 * pkin(5) + t204 * qJ(6) + t234) - g(3) * ((pkin(5) * t220 + qJ(6) * t219) * t221 + t237);];
U_reg  = t1;
