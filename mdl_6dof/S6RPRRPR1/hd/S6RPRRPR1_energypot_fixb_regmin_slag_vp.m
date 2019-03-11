% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:58:44
% EndTime: 2019-03-09 04:58:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (82->32), mult. (66->46), div. (0->0), fcn. (64->12), ass. (0->25)
t197 = qJ(3) + qJ(4);
t192 = pkin(11) + t197;
t211 = g(3) * sin(t192);
t210 = qJ(2) + pkin(6);
t196 = qJ(1) + pkin(10);
t190 = sin(t196);
t198 = sin(qJ(6));
t209 = t190 * t198;
t201 = cos(qJ(6));
t208 = t190 * t201;
t191 = cos(t196);
t207 = t191 * t198;
t206 = t191 * t201;
t205 = g(1) * t191 + g(2) * t190;
t200 = sin(qJ(1));
t203 = cos(qJ(1));
t204 = -g(1) * t203 - g(2) * t200;
t202 = cos(qJ(3));
t199 = sin(qJ(3));
t195 = -qJ(5) - pkin(8) - pkin(7);
t194 = cos(t197);
t193 = sin(t197);
t189 = cos(t192);
t187 = t202 * pkin(3) + pkin(4) * t194 + pkin(2);
t1 = [0, t204, g(1) * t200 - g(2) * t203, t204 * pkin(1) - g(3) * t210, 0, 0, 0, 0, 0, -g(3) * t199 - t205 * t202, -g(3) * t202 + t205 * t199, 0, 0, 0, 0, 0, -g(3) * t193 - t205 * t194, -g(3) * t194 + t205 * t193, -g(1) * t190 + g(2) * t191, -g(1) * (t203 * pkin(1) + t191 * t187 - t190 * t195) - g(2) * (t200 * pkin(1) + t190 * t187 + t191 * t195) - g(3) * (t199 * pkin(3) + pkin(4) * t193 + t210) 0, 0, 0, 0, 0, -g(1) * (t189 * t206 + t209) - g(2) * (t189 * t208 - t207) - t201 * t211, -g(1) * (-t189 * t207 + t208) - g(2) * (-t189 * t209 - t206) + t198 * t211;];
U_reg  = t1;
