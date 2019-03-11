% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:13:54
% EndTime: 2019-03-09 10:13:54
% DurationCPUTime: 0.07s
% Computational Cost: add. (104->39), mult. (91->48), div. (0->0), fcn. (89->10), ass. (0->25)
t207 = qJ(2) + pkin(10);
t202 = qJ(4) + t207;
t200 = cos(t202);
t222 = g(3) * t200;
t208 = -pkin(7) - qJ(3);
t210 = sin(qJ(2));
t221 = t210 * pkin(2) + pkin(6);
t213 = cos(qJ(2));
t201 = t213 * pkin(2) + pkin(1);
t209 = sin(qJ(6));
t211 = sin(qJ(1));
t220 = t211 * t209;
t212 = cos(qJ(6));
t219 = t211 * t212;
t214 = cos(qJ(1));
t218 = t214 * t209;
t217 = t214 * t212;
t216 = g(1) * t214 + g(2) * t211;
t199 = sin(t202);
t215 = pkin(4) * t200 + qJ(5) * t199 + pkin(3) * cos(t207) + t201;
t206 = -pkin(8) + t208;
t198 = g(1) * t211 - g(2) * t214;
t196 = g(3) * t199 + t216 * t200;
t195 = t216 * t199 - t222;
t1 = [0, -t216, t198, 0, 0, 0, 0, 0, -g(3) * t210 - t216 * t213, -g(3) * t213 + t216 * t210, -t198, -g(1) * (t214 * t201 - t211 * t208) - g(2) * (t211 * t201 + t214 * t208) - g(3) * t221, 0, 0, 0, 0, 0, -t196, t195, -t198, t196, -t195, -g(3) * (t199 * pkin(4) - t200 * qJ(5) + pkin(3) * sin(t207) + t221) + (-g(1) * t215 - g(2) * t206) * t214 + (g(1) * t206 - g(2) * t215) * t211, 0, 0, 0, 0, 0, -g(1) * (t199 * t218 + t219) - g(2) * (t199 * t220 - t217) + t209 * t222, -g(1) * (t199 * t217 - t220) - g(2) * (t199 * t219 + t218) + t212 * t222;];
U_reg  = t1;
