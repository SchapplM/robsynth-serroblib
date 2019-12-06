% Calculate inertial parameters regressor of potential energy for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:18
% EndTime: 2019-12-05 15:01:18
% DurationCPUTime: 0.18s
% Computational Cost: add. (131->61), mult. (143->79), div. (0->0), fcn. (138->10), ass. (0->36)
t58 = cos(pkin(9));
t45 = t58 * pkin(4) + pkin(3);
t54 = pkin(8) + qJ(3);
t48 = sin(t54);
t50 = cos(t54);
t61 = -pkin(6) - qJ(4);
t80 = t45 * t50 - t48 * t61;
t79 = g(3) * t48;
t78 = g(3) * qJ(1);
t53 = pkin(9) + qJ(5);
t47 = sin(t53);
t57 = sin(pkin(7));
t75 = t57 * t47;
t49 = cos(t53);
t74 = t57 * t49;
t55 = sin(pkin(9));
t73 = t57 * t55;
t72 = t57 * t58;
t60 = cos(pkin(7));
t71 = t60 * t47;
t70 = t60 * t49;
t69 = t60 * t55;
t68 = t60 * t58;
t59 = cos(pkin(8));
t46 = t59 * pkin(2) + pkin(1);
t62 = -pkin(5) - qJ(2);
t67 = t57 * t46 + t60 * t62;
t56 = sin(pkin(8));
t66 = t56 * pkin(2) + qJ(1);
t42 = t60 * t46;
t65 = -t57 * t62 + t42;
t64 = g(1) * t60 + g(2) * t57;
t63 = pkin(3) * t50 + qJ(4) * t48;
t40 = g(1) * t57 - g(2) * t60;
t39 = -g(3) * t50 + t64 * t48;
t1 = [0, 0, 0, 0, 0, 0, -t64, t40, -g(3), -t78, 0, 0, 0, 0, 0, 0, -g(3) * t56 - t64 * t59, -g(3) * t59 + t64 * t56, -t40, -g(1) * (t60 * pkin(1) + t57 * qJ(2)) - g(2) * (t57 * pkin(1) - t60 * qJ(2)) - t78, 0, 0, 0, 0, 0, 0, -t64 * t50 - t79, t39, -t40, -g(1) * t65 - g(2) * t67 - g(3) * t66, 0, 0, 0, 0, 0, 0, -g(1) * (t50 * t68 + t73) - g(2) * (t50 * t72 - t69) - t58 * t79, -g(1) * (-t50 * t69 + t72) - g(2) * (-t50 * t73 - t68) + t55 * t79, -t39, -g(1) * (t63 * t60 + t65) - g(2) * (t63 * t57 + t67) - g(3) * (t48 * pkin(3) - t50 * qJ(4) + t66), 0, 0, 0, 0, 0, 0, -g(1) * (t50 * t70 + t75) - g(2) * (t50 * t74 - t71) - t49 * t79, -g(1) * (-t50 * t71 + t74) - g(2) * (-t50 * t75 - t70) + t47 * t79, -t39, -g(1) * (t80 * t60 + t42) - g(2) * (-pkin(4) * t69 + t67) - g(3) * (t48 * t45 + t50 * t61 + t66) + (-g(1) * (pkin(4) * t55 - t62) - g(2) * t80) * t57;];
U_reg = t1;
