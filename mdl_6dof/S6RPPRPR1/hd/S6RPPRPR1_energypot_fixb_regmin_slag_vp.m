% Calculate minimal parameter regressor of potential energy for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:02
% EndTime: 2019-03-09 01:40:02
% DurationCPUTime: 0.11s
% Computational Cost: add. (134->49), mult. (107->70), div. (0->0), fcn. (106->12), ass. (0->33)
t192 = pkin(10) + qJ(4);
t184 = sin(t192);
t213 = g(3) * t184;
t198 = qJ(2) + pkin(6);
t212 = g(3) * t198;
t193 = qJ(1) + pkin(9);
t185 = sin(t193);
t187 = cos(t192);
t211 = t185 * t187;
t194 = sin(pkin(11));
t210 = t185 * t194;
t196 = cos(pkin(11));
t209 = t185 * t196;
t191 = pkin(11) + qJ(6);
t183 = sin(t191);
t188 = cos(t193);
t208 = t188 * t183;
t186 = cos(t191);
t207 = t188 * t186;
t206 = t188 * t194;
t205 = t188 * t196;
t204 = g(1) * t188 + g(2) * t185;
t200 = sin(qJ(1));
t201 = cos(qJ(1));
t203 = -g(1) * t201 - g(2) * t200;
t197 = cos(pkin(10));
t202 = t197 * pkin(3) + pkin(4) * t187 + qJ(5) * t184 + pkin(2);
t199 = -pkin(7) - qJ(3);
t195 = sin(pkin(10));
t190 = t201 * pkin(1);
t189 = t200 * pkin(1);
t181 = -g(3) * t187 + t204 * t184;
t1 = [0, t203, g(1) * t200 - g(2) * t201, pkin(1) * t203 - t212, -g(3) * t195 - t197 * t204, -g(3) * t197 + t195 * t204, -g(1) * t185 + g(2) * t188, -g(1) * (t188 * pkin(2) + t185 * qJ(3) + t190) - g(2) * (t185 * pkin(2) - t188 * qJ(3) + t189) - t212, 0, 0, 0, 0, 0, -t187 * t204 - t213, t181, -g(1) * (t187 * t205 + t210) - g(2) * (t187 * t209 - t206) - t196 * t213, -g(1) * (-t187 * t206 + t209) - g(2) * (-t187 * t210 - t205) + t194 * t213, -t181, -g(1) * t190 - g(2) * t189 - g(3) * (t195 * pkin(3) + t184 * pkin(4) - t187 * qJ(5) + t198) + (-g(1) * t202 - g(2) * t199) * t188 + (g(1) * t199 - g(2) * t202) * t185, 0, 0, 0, 0, 0, -g(1) * (t185 * t183 + t187 * t207) - g(2) * (t186 * t211 - t208) - t186 * t213, -g(1) * (t185 * t186 - t187 * t208) - g(2) * (-t183 * t211 - t207) + t183 * t213;];
U_reg  = t1;
