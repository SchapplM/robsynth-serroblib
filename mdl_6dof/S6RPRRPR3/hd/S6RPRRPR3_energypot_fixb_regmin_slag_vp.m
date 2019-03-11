% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:14
% EndTime: 2019-03-09 05:07:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (122->44), mult. (149->67), div. (0->0), fcn. (169->10), ass. (0->26)
t200 = sin(qJ(3));
t213 = g(3) * t200;
t212 = qJ(2) + pkin(6);
t199 = sin(qJ(4));
t204 = cos(qJ(3));
t211 = t199 * t204;
t203 = cos(qJ(4));
t210 = t203 * t204;
t197 = qJ(1) + pkin(10);
t195 = sin(t197);
t196 = cos(t197);
t209 = g(1) * t196 + g(2) * t195;
t201 = sin(qJ(1));
t205 = cos(qJ(1));
t208 = -g(1) * t205 - g(2) * t201;
t207 = pkin(3) * t204 + pkin(8) * t200 + pkin(2);
t190 = t195 * t211 + t196 * t203;
t192 = -t195 * t203 + t196 * t211;
t206 = g(1) * t192 + g(2) * t190 + t199 * t213;
t202 = cos(qJ(6));
t198 = sin(qJ(6));
t193 = t195 * t199 + t196 * t210;
t191 = t195 * t210 - t196 * t199;
t189 = -g(3) * t204 + t209 * t200;
t188 = -g(1) * t193 - g(2) * t191 - t203 * t213;
t1 = [0, t208, g(1) * t201 - g(2) * t205, t208 * pkin(1) - g(3) * t212, 0, 0, 0, 0, 0, -t209 * t204 - t213, t189, 0, 0, 0, 0, 0, t188, t206, t188, -t189, -t206, -g(1) * (t205 * pkin(1) + t193 * pkin(4) + t192 * qJ(5)) - g(2) * (t201 * pkin(1) + t191 * pkin(4) + t190 * qJ(5)) - g(3) * (-t204 * pkin(8) + t212) - (pkin(4) * t203 + qJ(5) * t199 + pkin(3)) * t213 + (g(2) * pkin(7) - g(1) * t207) * t196 + (-g(1) * pkin(7) - g(2) * t207) * t195, 0, 0, 0, 0, 0, -g(1) * (t192 * t198 + t193 * t202) - g(2) * (t190 * t198 + t191 * t202) - (t198 * t199 + t202 * t203) * t213, -g(1) * (t192 * t202 - t193 * t198) - g(2) * (t190 * t202 - t191 * t198) - (-t198 * t203 + t199 * t202) * t213;];
U_reg  = t1;
