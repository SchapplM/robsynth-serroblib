% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR2
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:02:33
% EndTime: 2019-03-09 05:02:33
% DurationCPUTime: 0.07s
% Computational Cost: add. (100->38), mult. (90->57), div. (0->0), fcn. (92->10), ass. (0->27)
t206 = sin(qJ(3));
t221 = g(3) * t206;
t220 = qJ(2) + pkin(6);
t203 = qJ(1) + pkin(10);
t200 = sin(t203);
t209 = cos(qJ(3));
t219 = t200 * t209;
t201 = cos(t203);
t218 = t201 * t209;
t205 = sin(qJ(4));
t217 = t205 * t209;
t208 = cos(qJ(4));
t216 = t208 * t209;
t215 = pkin(4) * t205 + pkin(7);
t214 = g(1) * t201 + g(2) * t200;
t207 = sin(qJ(1));
t210 = cos(qJ(1));
t213 = -g(1) * t210 - g(2) * t207;
t212 = t213 * pkin(1);
t199 = t208 * pkin(4) + pkin(3);
t204 = -qJ(5) - pkin(8);
t211 = t199 * t209 - t204 * t206 + pkin(2);
t202 = qJ(4) + pkin(11) + qJ(6);
t198 = cos(t202);
t197 = sin(t202);
t196 = -g(3) * t209 + t214 * t206;
t1 = [0, t213, g(1) * t207 - g(2) * t210, -g(3) * t220 + t212, 0, 0, 0, 0, 0, -t214 * t209 - t221, t196, 0, 0, 0, 0, 0, -g(1) * (t200 * t205 + t201 * t216) - g(2) * (t200 * t216 - t201 * t205) - t208 * t221, -g(1) * (t200 * t208 - t201 * t217) - g(2) * (-t200 * t217 - t201 * t208) + t205 * t221, -t196, -g(3) * (t206 * t199 + t209 * t204 + t220) + t212 + (-g(1) * t211 + g(2) * t215) * t201 + (-g(1) * t215 - g(2) * t211) * t200, 0, 0, 0, 0, 0, -g(1) * (t200 * t197 + t198 * t218) - g(2) * (-t201 * t197 + t198 * t219) - t198 * t221, -g(1) * (-t197 * t218 + t200 * t198) - g(2) * (-t197 * t219 - t201 * t198) + t197 * t221;];
U_reg  = t1;
