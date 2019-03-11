% Calculate minimal parameter regressor of potential energy for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:08
% EndTime: 2019-03-09 01:47:08
% DurationCPUTime: 0.09s
% Computational Cost: add. (74->37), mult. (120->54), div. (0->0), fcn. (137->10), ass. (0->27)
t202 = cos(qJ(1));
t201 = sin(qJ(1));
t183 = qJ(4) + pkin(10);
t200 = g(3) * sin(t183);
t199 = -qJ(3) + pkin(6);
t198 = cos(pkin(9));
t197 = sin(pkin(9));
t177 = cos(t183);
t185 = sin(qJ(6));
t196 = t177 * t185;
t187 = cos(qJ(6));
t195 = t177 * t187;
t194 = t202 * pkin(1) + t201 * qJ(2);
t193 = t202 * pkin(2) + t194;
t168 = -t201 * t197 - t202 * t198;
t169 = t202 * t197 - t201 * t198;
t192 = g(1) * t169 - g(2) * t168;
t191 = g(1) * t168 + g(2) * t169;
t190 = t201 * pkin(1) - t202 * qJ(2);
t189 = t201 * pkin(2) + t190;
t188 = cos(qJ(4));
t186 = sin(qJ(4));
t184 = -qJ(5) - pkin(7);
t175 = t188 * pkin(4) + pkin(3);
t171 = -g(1) * t202 - g(2) * t201;
t170 = g(1) * t201 - g(2) * t202;
t1 = [0, t171, t170, t171, -t170, -g(3) * pkin(6) - g(1) * t194 - g(2) * t190, t191, t192, -g(1) * t193 - g(2) * t189 - g(3) * t199, 0, 0, 0, 0, 0, g(3) * t186 + t191 * t188, g(3) * t188 - t191 * t186, -t192, -g(1) * (-t168 * t175 - t169 * t184 + t193) - g(2) * (t168 * t184 - t169 * t175 + t189) - g(3) * (-t186 * pkin(4) + t199) 0, 0, 0, 0, 0, -g(1) * (-t168 * t195 + t169 * t185) - g(2) * (-t168 * t185 - t169 * t195) + t187 * t200, -g(1) * (t168 * t196 + t169 * t187) - g(2) * (-t168 * t187 + t169 * t196) - t185 * t200;];
U_reg  = t1;
