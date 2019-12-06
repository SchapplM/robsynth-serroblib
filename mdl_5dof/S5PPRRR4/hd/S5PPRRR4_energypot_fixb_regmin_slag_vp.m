% Calculate minimal parameter regressor of potential energy for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:43
% EndTime: 2019-12-05 15:19:43
% DurationCPUTime: 0.10s
% Computational Cost: add. (136->47), mult. (373->93), div. (0->0), fcn. (487->14), ass. (0->40)
t203 = sin(pkin(10));
t209 = cos(pkin(5));
t224 = t203 * t209;
t204 = sin(pkin(6));
t205 = sin(pkin(5));
t223 = t204 * t205;
t222 = t204 * t209;
t207 = cos(pkin(10));
t221 = t205 * t207;
t208 = cos(pkin(6));
t220 = t205 * t208;
t206 = cos(pkin(11));
t219 = t206 * t208;
t218 = t207 * t209;
t202 = sin(pkin(11));
t198 = -t203 * t202 + t206 * t218;
t217 = -t198 * t208 + t204 * t221;
t200 = -t207 * t202 - t206 * t224;
t216 = t200 * t208 + t203 * t223;
t215 = cos(qJ(3));
t214 = cos(qJ(4));
t213 = cos(qJ(5));
t212 = sin(qJ(3));
t211 = sin(qJ(4));
t210 = sin(qJ(5));
t201 = -t202 * t224 + t207 * t206;
t199 = t202 * t218 + t203 * t206;
t197 = -t206 * t223 + t209 * t208;
t196 = -t200 * t204 + t203 * t220;
t195 = -t198 * t204 - t207 * t220;
t194 = t212 * t222 + (t202 * t215 + t212 * t219) * t205;
t193 = -t215 * t222 + (t202 * t212 - t215 * t219) * t205;
t192 = t194 * t214 + t197 * t211;
t191 = t201 * t215 + t216 * t212;
t190 = t201 * t212 - t216 * t215;
t189 = t199 * t215 - t217 * t212;
t188 = t199 * t212 + t217 * t215;
t187 = t191 * t214 + t196 * t211;
t186 = t189 * t214 + t195 * t211;
t1 = [-g(3) * qJ(1), -g(1) * (t203 * t205 * qJ(2) + t207 * pkin(1)) - g(2) * (t203 * pkin(1) - qJ(2) * t221) - g(3) * (t209 * qJ(2) + qJ(1)), 0, -g(1) * t191 - g(2) * t189 - g(3) * t194, g(1) * t190 + g(2) * t188 + g(3) * t193, 0, 0, 0, 0, 0, -g(1) * t187 - g(2) * t186 - g(3) * t192, -g(1) * (-t191 * t211 + t196 * t214) - g(2) * (-t189 * t211 + t195 * t214) - g(3) * (-t194 * t211 + t197 * t214), 0, 0, 0, 0, 0, -g(1) * (t187 * t213 + t190 * t210) - g(2) * (t186 * t213 + t188 * t210) - g(3) * (t192 * t213 + t193 * t210), -g(1) * (-t187 * t210 + t190 * t213) - g(2) * (-t186 * t210 + t188 * t213) - g(3) * (-t192 * t210 + t193 * t213);];
U_reg = t1;
