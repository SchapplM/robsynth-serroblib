% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP3
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
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:15
% EndTime: 2019-03-09 06:05:16
% DurationCPUTime: 0.09s
% Computational Cost: add. (144->46), mult. (130->67), div. (0->0), fcn. (140->10), ass. (0->32)
t200 = sin(qJ(3));
t216 = g(3) * t200;
t215 = qJ(2) + pkin(6);
t198 = qJ(4) + qJ(5);
t195 = sin(t198);
t203 = cos(qJ(3));
t214 = t195 * t203;
t196 = cos(t198);
t213 = t196 * t203;
t199 = sin(qJ(4));
t212 = t199 * t203;
t202 = cos(qJ(4));
t211 = t202 * t203;
t210 = pkin(4) * t199 + pkin(7);
t197 = qJ(1) + pkin(10);
t193 = sin(t197);
t194 = cos(t197);
t209 = g(1) * t194 + g(2) * t193;
t201 = sin(qJ(1));
t204 = cos(qJ(1));
t208 = -g(1) * t204 - g(2) * t201;
t192 = t202 * pkin(4) + pkin(3);
t205 = -pkin(9) - pkin(8);
t207 = t192 * t203 - t200 * t205 + pkin(2);
t186 = t193 * t214 + t194 * t196;
t188 = -t193 * t196 + t194 * t214;
t206 = g(1) * t188 + g(2) * t186 + t195 * t216;
t190 = -g(3) * t203 + t209 * t200;
t189 = t193 * t195 + t194 * t213;
t187 = t193 * t213 - t194 * t195;
t185 = -g(1) * t189 - g(2) * t187 - t196 * t216;
t1 = [0, t208, g(1) * t201 - g(2) * t204, t208 * pkin(1) - g(3) * t215, 0, 0, 0, 0, 0, -t209 * t203 - t216, t190, 0, 0, 0, 0, 0, -g(1) * (t193 * t199 + t194 * t211) - g(2) * (t193 * t211 - t194 * t199) - t202 * t216, -g(1) * (t193 * t202 - t194 * t212) - g(2) * (-t193 * t212 - t194 * t202) + t199 * t216, 0, 0, 0, 0, 0, t185, t206, t185, -t190, -t206, -g(1) * (t204 * pkin(1) + t189 * pkin(5) + t188 * qJ(6)) - g(2) * (t201 * pkin(1) + t187 * pkin(5) + t186 * qJ(6)) - g(3) * (t203 * t205 + t215) - (pkin(5) * t196 + qJ(6) * t195 + t192) * t216 + (-g(1) * t207 + g(2) * t210) * t194 + (-g(1) * t210 - g(2) * t207) * t193;];
U_reg  = t1;
