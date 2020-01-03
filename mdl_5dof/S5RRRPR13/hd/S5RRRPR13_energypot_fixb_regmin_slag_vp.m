% Calculate minimal parameter regressor of potential energy for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:32
% EndTime: 2019-12-31 21:46:32
% DurationCPUTime: 0.09s
% Computational Cost: add. (97->50), mult. (234->80), div. (0->0), fcn. (291->10), ass. (0->33)
t202 = sin(pkin(5));
t206 = sin(qJ(2));
t223 = t202 * t206;
t207 = sin(qJ(1));
t222 = t202 * t207;
t209 = cos(qJ(3));
t221 = t202 * t209;
t210 = cos(qJ(2));
t220 = t202 * t210;
t211 = cos(qJ(1));
t219 = t202 * t211;
t218 = t207 * t206;
t217 = t207 * t210;
t216 = t211 * t206;
t215 = t211 * t210;
t203 = cos(pkin(5));
t196 = t203 * t216 + t217;
t205 = sin(qJ(3));
t189 = t196 * t205 + t209 * t219;
t198 = -t203 * t218 + t215;
t191 = t198 * t205 - t207 * t221;
t193 = -t203 * t209 + t205 * t223;
t214 = g(1) * t191 + g(2) * t189 + g(3) * t193;
t190 = t196 * t209 - t205 * t219;
t192 = t198 * t209 + t205 * t222;
t194 = t203 * t205 + t206 * t221;
t213 = g(1) * t192 + g(2) * t190 + g(3) * t194;
t195 = -t203 * t215 + t218;
t197 = t203 * t217 + t216;
t212 = -g(1) * t197 - g(2) * t195 + g(3) * t220;
t208 = cos(qJ(5));
t204 = sin(qJ(5));
t1 = [0, -g(1) * t211 - g(2) * t207, g(1) * t207 - g(2) * t211, 0, 0, 0, 0, 0, -g(1) * t198 - g(2) * t196 - g(3) * t223, -t212, 0, 0, 0, 0, 0, -t213, t214, t212, t213, -t214, -g(1) * (t211 * pkin(1) + t198 * pkin(2) + t192 * pkin(3) + pkin(7) * t222 + t197 * pkin(8) + t191 * qJ(4)) - g(2) * (t207 * pkin(1) + t196 * pkin(2) + t190 * pkin(3) - pkin(7) * t219 + t195 * pkin(8) + t189 * qJ(4)) - g(3) * (t194 * pkin(3) + t203 * pkin(7) + t193 * qJ(4) + pkin(6) + (pkin(2) * t206 - pkin(8) * t210) * t202), 0, 0, 0, 0, 0, -g(1) * (t191 * t204 + t197 * t208) - g(2) * (t189 * t204 + t195 * t208) - g(3) * (t193 * t204 - t208 * t220), -g(1) * (t191 * t208 - t197 * t204) - g(2) * (t189 * t208 - t195 * t204) - g(3) * (t193 * t208 + t204 * t220);];
U_reg = t1;
