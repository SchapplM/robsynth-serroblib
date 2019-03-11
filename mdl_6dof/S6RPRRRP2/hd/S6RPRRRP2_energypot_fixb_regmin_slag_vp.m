% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP2
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
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:15
% EndTime: 2019-03-09 06:01:15
% DurationCPUTime: 0.09s
% Computational Cost: add. (103->40), mult. (95->59), div. (0->0), fcn. (97->10), ass. (0->27)
t190 = sin(qJ(3));
t205 = g(3) * t190;
t188 = qJ(4) + qJ(5);
t184 = sin(t188);
t189 = sin(qJ(4));
t204 = pkin(4) * t189 + pkin(5) * t184 + pkin(7);
t203 = qJ(2) + pkin(6);
t193 = cos(qJ(3));
t202 = t184 * t193;
t185 = cos(t188);
t201 = t185 * t193;
t200 = t189 * t193;
t192 = cos(qJ(4));
t199 = t192 * t193;
t187 = qJ(1) + pkin(10);
t182 = sin(t187);
t183 = cos(t187);
t198 = g(1) * t183 + g(2) * t182;
t191 = sin(qJ(1));
t194 = cos(qJ(1));
t197 = -g(1) * t194 - g(2) * t191;
t196 = t197 * pkin(1);
t180 = pkin(4) * t192 + pkin(5) * t185 + pkin(3);
t186 = -qJ(6) - pkin(9) - pkin(8);
t195 = t180 * t193 - t186 * t190 + pkin(2);
t179 = -g(3) * t193 + t190 * t198;
t1 = [0, t197, g(1) * t191 - g(2) * t194, -g(3) * t203 + t196, 0, 0, 0, 0, 0, -t193 * t198 - t205, t179, 0, 0, 0, 0, 0, -g(1) * (t182 * t189 + t183 * t199) - g(2) * (t182 * t199 - t183 * t189) - t192 * t205, -g(1) * (t182 * t192 - t183 * t200) - g(2) * (-t182 * t200 - t183 * t192) + t189 * t205, 0, 0, 0, 0, 0, -g(1) * (t182 * t184 + t183 * t201) - g(2) * (t182 * t201 - t183 * t184) - t185 * t205, -g(1) * (t182 * t185 - t183 * t202) - g(2) * (-t182 * t202 - t183 * t185) + t184 * t205, -t179, -g(3) * (t180 * t190 + t186 * t193 + t203) + t196 + (-g(1) * t195 + g(2) * t204) * t183 + (-g(1) * t204 - g(2) * t195) * t182;];
U_reg  = t1;
