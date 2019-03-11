% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:26
% EndTime: 2019-03-09 03:42:26
% DurationCPUTime: 0.08s
% Computational Cost: add. (116->44), mult. (103->69), div. (0->0), fcn. (109->12), ass. (0->27)
t209 = sin(qJ(3));
t222 = g(3) * t209;
t221 = qJ(2) + pkin(6);
t206 = qJ(1) + pkin(10);
t201 = sin(t206);
t211 = cos(qJ(3));
t220 = t201 * t211;
t203 = cos(t206);
t219 = t203 * t211;
t207 = sin(pkin(11));
t218 = t207 * t211;
t208 = cos(pkin(11));
t217 = t208 * t211;
t205 = pkin(11) + qJ(5);
t216 = g(1) * t203 + g(2) * t201;
t210 = sin(qJ(1));
t212 = cos(qJ(1));
t215 = -g(1) * t212 - g(2) * t210;
t214 = pkin(3) * t211 + qJ(4) * t209 + pkin(2);
t213 = t215 * pkin(1);
t204 = qJ(6) + t205;
t202 = cos(t205);
t200 = sin(t205);
t199 = cos(t204);
t198 = sin(t204);
t197 = -g(3) * t211 + t216 * t209;
t1 = [0, t215, g(1) * t210 - g(2) * t212, -g(3) * t221 + t213, 0, 0, 0, 0, 0, -t216 * t211 - t222, t197, -g(1) * (t201 * t207 + t203 * t217) - g(2) * (t201 * t217 - t203 * t207) - t208 * t222, -g(1) * (t201 * t208 - t203 * t218) - g(2) * (-t201 * t218 - t203 * t208) + t207 * t222, -t197, -g(3) * (t209 * pkin(3) - t211 * qJ(4) + t221) + t213 + (g(2) * pkin(7) - g(1) * t214) * t203 + (-g(1) * pkin(7) - g(2) * t214) * t201, 0, 0, 0, 0, 0, -g(1) * (t201 * t200 + t202 * t219) - g(2) * (-t203 * t200 + t202 * t220) - t202 * t222, -g(1) * (-t200 * t219 + t201 * t202) - g(2) * (-t200 * t220 - t203 * t202) + t200 * t222, 0, 0, 0, 0, 0, -g(1) * (t201 * t198 + t199 * t219) - g(2) * (-t203 * t198 + t199 * t220) - t199 * t222, -g(1) * (-t198 * t219 + t201 * t199) - g(2) * (-t198 * t220 - t203 * t199) + t198 * t222;];
U_reg  = t1;
