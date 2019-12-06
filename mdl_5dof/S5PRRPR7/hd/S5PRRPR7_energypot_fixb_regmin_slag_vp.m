% Calculate minimal parameter regressor of potential energy for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:33
% EndTime: 2019-12-05 16:37:33
% DurationCPUTime: 0.15s
% Computational Cost: add. (125->57), mult. (307->97), div. (0->0), fcn. (392->12), ass. (0->34)
t190 = sin(pkin(5));
t207 = pkin(6) * t190;
t195 = sin(qJ(3));
t206 = t190 * t195;
t196 = sin(qJ(2));
t205 = t190 * t196;
t198 = cos(qJ(3));
t204 = t190 * t198;
t199 = cos(qJ(2));
t203 = t190 * t199;
t193 = cos(pkin(5));
t202 = t193 * t196;
t201 = t193 * t199;
t189 = sin(pkin(9));
t192 = cos(pkin(9));
t180 = t189 * t199 + t192 * t202;
t175 = t180 * t195 + t192 * t204;
t182 = -t189 * t202 + t192 * t199;
t177 = t182 * t195 - t189 * t204;
t183 = -t193 * t198 + t195 * t205;
t200 = g(1) * t177 + g(2) * t175 + g(3) * t183;
t197 = cos(qJ(5));
t194 = sin(qJ(5));
t191 = cos(pkin(10));
t188 = sin(pkin(10));
t184 = t193 * t195 + t196 * t204;
t181 = t189 * t201 + t192 * t196;
t179 = t189 * t196 - t192 * t201;
t178 = t182 * t198 + t189 * t206;
t176 = t180 * t198 - t192 * t206;
t174 = t184 * t191 - t188 * t203;
t173 = t178 * t191 + t181 * t188;
t172 = t176 * t191 + t179 * t188;
t1 = [-g(3) * qJ(1), 0, -g(1) * t182 - g(2) * t180 - g(3) * t205, g(1) * t181 + g(2) * t179 - g(3) * t203, 0, 0, 0, 0, 0, -g(1) * t178 - g(2) * t176 - g(3) * t184, t200, -g(1) * t173 - g(2) * t172 - g(3) * t174, -g(1) * (-t178 * t188 + t181 * t191) - g(2) * (-t176 * t188 + t179 * t191) - g(3) * (-t184 * t188 - t191 * t203), -t200, -g(1) * (t192 * pkin(1) + t182 * pkin(2) + t178 * pkin(3) + t181 * pkin(7) + t177 * qJ(4) + t189 * t207) - g(2) * (t189 * pkin(1) + t180 * pkin(2) + t176 * pkin(3) + t179 * pkin(7) + t175 * qJ(4) - t192 * t207) - g(3) * (t184 * pkin(3) + t193 * pkin(6) + t183 * qJ(4) + qJ(1) + (pkin(2) * t196 - pkin(7) * t199) * t190), 0, 0, 0, 0, 0, -g(1) * (t173 * t197 + t177 * t194) - g(2) * (t172 * t197 + t175 * t194) - g(3) * (t174 * t197 + t183 * t194), -g(1) * (-t173 * t194 + t177 * t197) - g(2) * (-t172 * t194 + t175 * t197) - g(3) * (-t174 * t194 + t183 * t197);];
U_reg = t1;
