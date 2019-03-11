% Calculate minimal parameter regressor of potential energy for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:05
% EndTime: 2019-03-08 21:21:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (183->70), mult. (405->105), div. (0->0), fcn. (504->12), ass. (0->38)
t251 = pkin(4) + pkin(8);
t230 = sin(pkin(6));
t250 = pkin(7) * t230;
t234 = sin(qJ(3));
t249 = t230 * t234;
t235 = sin(qJ(2));
t248 = t230 * t235;
t236 = cos(qJ(3));
t247 = t230 * t236;
t237 = cos(qJ(2));
t246 = t230 * t237;
t233 = cos(pkin(6));
t245 = t233 * t235;
t244 = t233 * t237;
t214 = -t233 * t236 + t234 * t248;
t215 = t233 * t234 + t235 * t247;
t243 = pkin(2) * t248 + t215 * pkin(3) + t233 * pkin(7) + t214 * qJ(4) + qJ(1);
t229 = sin(pkin(10));
t232 = cos(pkin(10));
t213 = -t229 * t245 + t232 * t237;
t205 = t213 * t234 - t229 * t247;
t206 = t213 * t236 + t229 * t249;
t242 = t232 * pkin(1) + t213 * pkin(2) + t206 * pkin(3) + t205 * qJ(4) + t229 * t250;
t211 = t229 * t237 + t232 * t245;
t203 = t211 * t234 + t232 * t247;
t241 = g(1) * t205 + g(2) * t203 + g(3) * t214;
t204 = t211 * t236 - t232 * t249;
t240 = g(1) * t206 + g(2) * t204 + g(3) * t215;
t210 = t229 * t235 - t232 * t244;
t212 = t229 * t244 + t232 * t235;
t239 = -g(1) * t212 - g(2) * t210 + g(3) * t246;
t238 = t229 * pkin(1) + t211 * pkin(2) + t204 * pkin(3) + t203 * qJ(4) - t232 * t250;
t231 = cos(pkin(11));
t228 = sin(pkin(11));
t227 = pkin(11) + qJ(6);
t223 = cos(t227);
t222 = sin(t227);
t1 = [-g(3) * qJ(1), 0, -g(1) * t213 - g(2) * t211 - g(3) * t248, -t239, 0, 0, 0, 0, 0, -t240, t241, t239, t240, -t241, -g(1) * (t212 * pkin(8) + t242) - g(2) * (t210 * pkin(8) + t238) - g(3) * (-pkin(8) * t246 + t243) -g(1) * (t205 * t228 + t212 * t231) - g(2) * (t203 * t228 + t210 * t231) - g(3) * (t214 * t228 - t231 * t246) -g(1) * (t205 * t231 - t212 * t228) - g(2) * (t203 * t231 - t210 * t228) - g(3) * (t214 * t231 + t228 * t246) -t240, -g(1) * (t206 * qJ(5) + t251 * t212 + t242) - g(2) * (t204 * qJ(5) + t251 * t210 + t238) - g(3) * (t215 * qJ(5) - t251 * t246 + t243) 0, 0, 0, 0, 0, -g(1) * (t205 * t222 + t212 * t223) - g(2) * (t203 * t222 + t210 * t223) - g(3) * (t214 * t222 - t223 * t246) -g(1) * (t205 * t223 - t212 * t222) - g(2) * (t203 * t223 - t210 * t222) - g(3) * (t214 * t223 + t222 * t246);];
U_reg  = t1;
