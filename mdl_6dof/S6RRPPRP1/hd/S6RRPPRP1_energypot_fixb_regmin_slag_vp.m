% Calculate minimal parameter regressor of potential energy for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:47
% EndTime: 2019-03-09 08:27:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (164->57), mult. (163->76), div. (0->0), fcn. (170->10), ass. (0->40)
t222 = qJ(2) + pkin(9);
t217 = sin(t222);
t246 = g(3) * t217;
t227 = sin(qJ(2));
t245 = t227 * pkin(2) + pkin(6);
t224 = cos(pkin(10));
t213 = t224 * pkin(4) + pkin(3);
t219 = cos(t222);
t244 = t213 * t219;
t221 = pkin(10) + qJ(5);
t216 = sin(t221);
t228 = sin(qJ(1));
t243 = t228 * t216;
t218 = cos(t221);
t242 = t228 * t218;
t223 = sin(pkin(10));
t241 = t228 * t223;
t240 = t228 * t224;
t230 = cos(qJ(1));
t239 = t230 * t216;
t238 = t230 * t218;
t237 = t230 * t223;
t236 = t230 * t224;
t229 = cos(qJ(2));
t215 = t229 * pkin(2) + pkin(1);
t226 = -pkin(7) - qJ(3);
t235 = t228 * t215 + t230 * t226;
t234 = t230 * t215 - t228 * t226;
t233 = g(1) * t230 + g(2) * t228;
t232 = pkin(3) * t219 + qJ(4) * t217;
t205 = t219 * t243 + t238;
t207 = t219 * t239 - t242;
t231 = g(1) * t207 + g(2) * t205 + t216 * t246;
t225 = -pkin(8) - qJ(4);
t209 = g(1) * t228 - g(2) * t230;
t208 = t219 * t238 + t243;
t206 = t219 * t242 - t239;
t204 = g(3) * t219 - t233 * t217;
t203 = -g(1) * t208 - g(2) * t206 - t218 * t246;
t1 = [0, -t233, t209, 0, 0, 0, 0, 0, -g(3) * t227 - t233 * t229, -g(3) * t229 + t233 * t227, -t209, -g(1) * t234 - g(2) * t235 - g(3) * t245, -g(1) * (t219 * t236 + t241) - g(2) * (t219 * t240 - t237) - t224 * t246, -g(1) * (-t219 * t237 + t240) - g(2) * (-t219 * t241 - t236) + t223 * t246, t204, -g(1) * (t232 * t230 + t234) - g(2) * (t232 * t228 + t235) - g(3) * (t217 * pkin(3) - t219 * qJ(4) + t245) 0, 0, 0, 0, 0, t203, t231, t203, t204, -t231, -g(1) * (pkin(4) * t241 + t208 * pkin(5) + t207 * qJ(6) + t230 * t244 + t234) - g(2) * (-pkin(4) * t237 + t206 * pkin(5) + t205 * qJ(6) + t228 * t244 + t235) - g(3) * (t219 * t225 + t245) + (-g(3) * (pkin(5) * t218 + qJ(6) * t216 + t213) + t233 * t225) * t217;];
U_reg  = t1;
