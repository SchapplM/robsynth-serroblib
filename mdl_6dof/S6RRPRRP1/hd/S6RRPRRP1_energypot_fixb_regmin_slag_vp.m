% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:14
% EndTime: 2019-03-09 11:41:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (103->40), mult. (91->50), div. (0->0), fcn. (89->10), ass. (0->26)
t205 = qJ(2) + pkin(10);
t201 = qJ(4) + t205;
t197 = sin(t201);
t222 = g(3) * t197;
t207 = -pkin(7) - qJ(3);
t209 = sin(qJ(2));
t221 = t209 * pkin(2) + pkin(6);
t212 = cos(qJ(2));
t200 = t212 * pkin(2) + pkin(1);
t208 = sin(qJ(5));
t210 = sin(qJ(1));
t220 = t210 * t208;
t211 = cos(qJ(5));
t219 = t210 * t211;
t213 = cos(qJ(1));
t218 = t213 * t208;
t217 = t213 * t211;
t216 = pkin(5) * t208 + pkin(8) - t207;
t215 = g(1) * t213 + g(2) * t210;
t198 = cos(t201);
t199 = t211 * pkin(5) + pkin(4);
t206 = -qJ(6) - pkin(9);
t214 = t197 * t206 - t198 * t199 - pkin(3) * cos(t205) - t200;
t196 = g(1) * t210 - g(2) * t213;
t194 = -g(3) * t198 + t215 * t197;
t1 = [0, -t215, t196, 0, 0, 0, 0, 0, -g(3) * t209 - t215 * t212, -g(3) * t212 + t215 * t209, -t196, -g(1) * (t213 * t200 - t210 * t207) - g(2) * (t210 * t200 + t213 * t207) - g(3) * t221, 0, 0, 0, 0, 0, -t215 * t198 - t222, t194, 0, 0, 0, 0, 0, -g(1) * (t198 * t217 + t220) - g(2) * (t198 * t219 - t218) - t211 * t222, -g(1) * (-t198 * t218 + t219) - g(2) * (-t198 * t220 - t217) + t208 * t222, -t194, -g(3) * (t197 * t199 + t198 * t206 + pkin(3) * sin(t205) + t221) + (g(1) * t214 + g(2) * t216) * t213 + (-g(1) * t216 + g(2) * t214) * t210;];
U_reg  = t1;
