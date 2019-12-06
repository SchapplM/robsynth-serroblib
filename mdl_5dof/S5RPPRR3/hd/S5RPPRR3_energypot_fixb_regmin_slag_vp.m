% Calculate minimal parameter regressor of potential energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:42:04
% EndTime: 2019-12-05 17:42:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (61->21), mult. (50->30), div. (0->0), fcn. (44->10), ass. (0->17)
t144 = g(1) * (qJ(2) + pkin(5));
t135 = pkin(9) + qJ(4);
t136 = qJ(1) + pkin(8);
t131 = sin(t136);
t133 = cos(t136);
t143 = g(2) * t131 - g(3) * t133;
t140 = sin(qJ(1));
t141 = cos(qJ(1));
t142 = g(2) * t140 - g(3) * t141;
t138 = cos(pkin(9));
t137 = sin(pkin(9));
t134 = qJ(5) + t135;
t132 = cos(t135);
t130 = sin(t135);
t129 = cos(t134);
t128 = sin(t134);
t1 = [0, t142, g(2) * t141 + g(3) * t140, t142 * pkin(1) - t144, -g(1) * t137 + t143 * t138, -g(1) * t138 - t143 * t137, -g(2) * t133 - g(3) * t131, -t144 - g(2) * (-t140 * pkin(1) - t131 * pkin(2) + t133 * qJ(3)) - g(3) * (t141 * pkin(1) + t133 * pkin(2) + t131 * qJ(3)), 0, 0, 0, 0, 0, -g(1) * t130 + t143 * t132, -g(1) * t132 - t143 * t130, 0, 0, 0, 0, 0, -g(1) * t128 + t143 * t129, -g(1) * t129 - t143 * t128;];
U_reg = t1;
