% Calculate minimal parameter regressor of potential energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:31
% EndTime: 2019-12-05 17:31:31
% DurationCPUTime: 0.12s
% Computational Cost: add. (93->57), mult. (210->89), div. (0->0), fcn. (241->10), ass. (0->37)
t195 = cos(pkin(7));
t215 = pkin(2) * t195;
t197 = sin(qJ(1));
t214 = g(2) * t197;
t191 = sin(pkin(8));
t192 = sin(pkin(7));
t213 = t191 * t192;
t194 = cos(pkin(8));
t212 = t192 * t194;
t211 = t192 * t197;
t199 = cos(qJ(1));
t210 = t192 * t199;
t209 = t197 * t191;
t208 = t197 * t194;
t207 = t199 * t191;
t206 = t199 * t194;
t205 = t199 * pkin(1) + t197 * qJ(2);
t204 = qJ(3) * t210 + t199 * t215 + t205;
t203 = t192 * pkin(2) - t195 * qJ(3) + pkin(5);
t202 = -g(3) * t199 + t214;
t201 = (-qJ(3) * t192 - pkin(1) - t215) * t214;
t178 = t195 * t209 + t206;
t180 = t195 * t207 - t208;
t200 = g(1) * t213 - g(2) * t178 + g(3) * t180;
t198 = cos(qJ(5));
t196 = sin(qJ(5));
t193 = cos(pkin(9));
t190 = sin(pkin(9));
t188 = t199 * qJ(2);
t182 = g(2) * t199 + g(3) * t197;
t181 = t195 * t206 + t209;
t179 = -t195 * t208 + t207;
t177 = -t195 * t190 + t193 * t212;
t176 = g(1) * t195 + t202 * t192;
t175 = t181 * t193 + t190 * t210;
t174 = t179 * t193 - t190 * t211;
t1 = [0, t202, t182, -g(1) * t192 + t202 * t195, -t176, -t182, -g(1) * pkin(5) - g(2) * (-t197 * pkin(1) + t188) - g(3) * t205, -g(1) * t212 - g(2) * t179 - g(3) * t181, t200, t176, -g(1) * t203 - g(2) * t188 - g(3) * t204 - t201, -g(1) * t177 - g(2) * t174 - g(3) * t175, -g(1) * (-t190 * t212 - t195 * t193) - g(2) * (-t179 * t190 - t193 * t211) - g(3) * (-t181 * t190 + t193 * t210), -t200, -g(1) * ((pkin(3) * t194 + qJ(4) * t191) * t192 + t203) - g(2) * (t179 * pkin(3) - t178 * qJ(4) + t188) - g(3) * (t181 * pkin(3) + t180 * qJ(4) + t204) - t201, 0, 0, 0, 0, 0, -g(1) * (t177 * t198 + t196 * t213) - g(2) * (t174 * t198 - t178 * t196) - g(3) * (t175 * t198 + t180 * t196), -g(1) * (-t177 * t196 + t198 * t213) - g(2) * (-t174 * t196 - t178 * t198) - g(3) * (-t175 * t196 + t180 * t198);];
U_reg = t1;
