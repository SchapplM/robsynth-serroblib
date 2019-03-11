% Calculate minimal parameter regressor of potential energy for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:36:46
% EndTime: 2019-03-08 23:36:46
% DurationCPUTime: 0.14s
% Computational Cost: add. (185->64), mult. (457->103), div. (0->0), fcn. (592->12), ass. (0->39)
t234 = sin(pkin(6));
t253 = pkin(7) * t234;
t239 = sin(qJ(3));
t252 = t234 * t239;
t240 = sin(qJ(2));
t251 = t234 * t240;
t243 = cos(qJ(3));
t250 = t234 * t243;
t244 = cos(qJ(2));
t249 = t234 * t244;
t236 = cos(pkin(6));
t248 = t236 * t240;
t247 = t236 * t244;
t233 = sin(pkin(11));
t235 = cos(pkin(11));
t225 = t233 * t244 + t235 * t248;
t217 = t225 * t243 - t235 * t252;
t224 = t233 * t240 - t235 * t247;
t238 = sin(qJ(4));
t242 = cos(qJ(4));
t212 = t217 * t238 - t224 * t242;
t227 = -t233 * t248 + t235 * t244;
t219 = t227 * t243 + t233 * t252;
t226 = t233 * t247 + t235 * t240;
t214 = t219 * t238 - t226 * t242;
t229 = t236 * t239 + t240 * t250;
t220 = t229 * t238 + t242 * t249;
t246 = g(1) * t214 + g(2) * t212 + g(3) * t220;
t216 = t225 * t239 + t235 * t250;
t218 = t227 * t239 - t233 * t250;
t228 = -t236 * t243 + t239 * t251;
t245 = g(1) * t218 + g(2) * t216 + g(3) * t228;
t241 = cos(qJ(6));
t237 = sin(qJ(6));
t221 = t229 * t242 - t238 * t249;
t215 = t219 * t242 + t226 * t238;
t213 = t217 * t242 + t224 * t238;
t211 = -g(1) * t215 - g(2) * t213 - g(3) * t221;
t1 = [-g(3) * qJ(1), 0, -g(1) * t227 - g(2) * t225 - g(3) * t251, g(1) * t226 + g(2) * t224 - g(3) * t249, 0, 0, 0, 0, 0, -g(1) * t219 - g(2) * t217 - g(3) * t229, t245, 0, 0, 0, 0, 0, t211, t246, t211, -t245, -t246, -g(1) * (t235 * pkin(1) + t227 * pkin(2) + t219 * pkin(3) + t215 * pkin(4) + t226 * pkin(8) + t218 * pkin(9) + t214 * qJ(5) + t233 * t253) - g(2) * (t233 * pkin(1) + t225 * pkin(2) + t217 * pkin(3) + t213 * pkin(4) + t224 * pkin(8) + t216 * pkin(9) + t212 * qJ(5) - t235 * t253) - g(3) * (t229 * pkin(3) + t221 * pkin(4) + t236 * pkin(7) + t228 * pkin(9) + t220 * qJ(5) + qJ(1) + (pkin(2) * t240 - pkin(8) * t244) * t234) 0, 0, 0, 0, 0, -g(1) * (t214 * t237 + t215 * t241) - g(2) * (t212 * t237 + t213 * t241) - g(3) * (t220 * t237 + t221 * t241) -g(1) * (t214 * t241 - t215 * t237) - g(2) * (t212 * t241 - t213 * t237) - g(3) * (t220 * t241 - t221 * t237);];
U_reg  = t1;
