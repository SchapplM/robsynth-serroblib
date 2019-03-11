% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP1
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:05
% EndTime: 2019-03-09 05:57:06
% DurationCPUTime: 0.09s
% Computational Cost: add. (134->40), mult. (116->55), div. (0->0), fcn. (122->10), ass. (0->28)
t188 = qJ(3) + qJ(4);
t185 = sin(t188);
t203 = g(3) * t185;
t202 = qJ(2) + pkin(6);
t186 = cos(t188);
t189 = sin(qJ(5));
t201 = t186 * t189;
t192 = cos(qJ(5));
t200 = t186 * t192;
t187 = qJ(1) + pkin(10);
t183 = sin(t187);
t184 = cos(t187);
t199 = g(1) * t184 + g(2) * t183;
t191 = sin(qJ(1));
t194 = cos(qJ(1));
t198 = -g(1) * t194 - g(2) * t191;
t193 = cos(qJ(3));
t197 = t193 * pkin(3) + pkin(4) * t186 + pkin(9) * t185 + pkin(2);
t177 = t183 * t201 + t184 * t192;
t179 = -t183 * t192 + t184 * t201;
t196 = g(1) * t179 + g(2) * t177 + t189 * t203;
t195 = -pkin(8) - pkin(7);
t190 = sin(qJ(3));
t180 = t183 * t189 + t184 * t200;
t178 = t183 * t200 - t184 * t189;
t176 = -g(3) * t186 + t199 * t185;
t175 = -g(1) * t180 - g(2) * t178 - t192 * t203;
t1 = [0, t198, g(1) * t191 - g(2) * t194, t198 * pkin(1) - g(3) * t202, 0, 0, 0, 0, 0, -g(3) * t190 - t199 * t193, -g(3) * t193 + t199 * t190, 0, 0, 0, 0, 0, -t199 * t186 - t203, t176, 0, 0, 0, 0, 0, t175, t196, t175, -t176, -t196, -g(1) * (t194 * pkin(1) + t180 * pkin(5) + t179 * qJ(6)) - g(2) * (t191 * pkin(1) + t178 * pkin(5) + t177 * qJ(6)) - g(3) * (t190 * pkin(3) - t186 * pkin(9) + t202) - (pkin(5) * t192 + qJ(6) * t189 + pkin(4)) * t203 + (-g(1) * t197 - g(2) * t195) * t184 + (g(1) * t195 - g(2) * t197) * t183;];
U_reg  = t1;
