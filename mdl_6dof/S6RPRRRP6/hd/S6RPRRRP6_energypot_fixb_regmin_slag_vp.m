% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:16:43
% EndTime: 2019-03-09 06:16:44
% DurationCPUTime: 0.10s
% Computational Cost: add. (106->47), mult. (110->64), div. (0->0), fcn. (112->10), ass. (0->29)
t201 = pkin(10) + qJ(3);
t196 = sin(t201);
t221 = g(3) * t196;
t202 = qJ(4) + qJ(5);
t198 = sin(t202);
t209 = cos(qJ(1));
t220 = t198 * t209;
t199 = cos(t202);
t219 = t199 * t209;
t206 = sin(qJ(4));
t218 = t206 * t209;
t207 = sin(qJ(1));
t217 = t207 * t198;
t216 = t207 * t199;
t215 = t207 * t206;
t208 = cos(qJ(4));
t214 = t207 * t208;
t213 = t208 * t209;
t212 = pkin(4) * t206 + pkin(5) * t198 + pkin(7) + qJ(2);
t211 = g(1) * t209 + g(2) * t207;
t192 = pkin(4) * t208 + pkin(5) * t199 + pkin(3);
t197 = cos(t201);
t200 = -qJ(6) - pkin(9) - pkin(8);
t204 = cos(pkin(10));
t210 = pkin(2) * t204 + t192 * t197 - t196 * t200 + pkin(1);
t203 = sin(pkin(10));
t194 = g(1) * t207 - g(2) * t209;
t191 = -g(3) * t197 + t211 * t196;
t1 = [0, -t211, t194, -g(3) * t203 - t211 * t204, -g(3) * t204 + t211 * t203, -t194, -g(1) * (pkin(1) * t209 + t207 * qJ(2)) - g(2) * (t207 * pkin(1) - qJ(2) * t209) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t211 * t197 - t221, t191, 0, 0, 0, 0, 0, -g(1) * (t197 * t213 + t215) - g(2) * (t197 * t214 - t218) - t208 * t221, -g(1) * (-t197 * t218 + t214) - g(2) * (-t197 * t215 - t213) + t206 * t221, 0, 0, 0, 0, 0, -g(1) * (t197 * t219 + t217) - g(2) * (t197 * t216 - t220) - t199 * t221, -g(1) * (-t197 * t220 + t216) - g(2) * (-t197 * t217 - t219) + t198 * t221, -t191, -g(3) * (pkin(2) * t203 + t192 * t196 + t197 * t200 + pkin(6)) + (-g(1) * t210 + g(2) * t212) * t209 + (-g(1) * t212 - g(2) * t210) * t207;];
U_reg  = t1;
