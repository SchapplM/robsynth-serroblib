% Calculate inertial parameters regressor of potential energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:29:04
% EndTime: 2019-12-05 17:29:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (138->58), mult. (130->79), div. (0->0), fcn. (125->10), ass. (0->32)
t59 = sin(pkin(8));
t79 = g(1) * t59;
t62 = qJ(2) + pkin(5);
t78 = g(1) * t62;
t57 = qJ(1) + pkin(7);
t52 = sin(t57);
t77 = g(2) * t52;
t61 = cos(pkin(8));
t76 = t52 * t61;
t54 = cos(t57);
t58 = sin(pkin(9));
t75 = t54 * t58;
t74 = t54 * t61;
t73 = t58 * t61;
t63 = -pkin(6) - qJ(4);
t72 = t59 * t63;
t60 = cos(pkin(9));
t71 = t60 * t61;
t65 = cos(qJ(1));
t70 = t65 * pkin(1) + t54 * pkin(2) + t52 * qJ(3);
t64 = sin(qJ(1));
t69 = -t64 * pkin(1) + t54 * qJ(3);
t68 = -g(3) * t54 + t77;
t67 = g(2) * t64 - g(3) * t65;
t66 = pkin(3) * t61 + qJ(4) * t59;
t56 = pkin(9) + qJ(5);
t53 = cos(t56);
t51 = sin(t56);
t50 = t60 * pkin(4) + pkin(3);
t46 = g(2) * t54 + g(3) * t52;
t45 = g(1) * t61 + t68 * t59;
t1 = [0, 0, 0, 0, 0, 0, t67, g(2) * t65 + g(3) * t64, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, t68, t46, -g(1), t67 * pkin(1) - t78, 0, 0, 0, 0, 0, 0, t68 * t61 - t79, -t45, -t46, -t78 - g(2) * (-t52 * pkin(2) + t69) - g(3) * t70, 0, 0, 0, 0, 0, 0, -t60 * t79 - g(2) * (-t52 * t71 + t75) - g(3) * (t52 * t58 + t54 * t71), t58 * t79 - g(2) * (t52 * t73 + t54 * t60) - g(3) * (t52 * t60 - t54 * t73), t45, -g(1) * (t59 * pkin(3) - t61 * qJ(4) + t62) - g(2) * t69 - g(3) * (t66 * t54 + t70) - (-pkin(2) - t66) * t77, 0, 0, 0, 0, 0, 0, -t53 * t79 - g(2) * (t54 * t51 - t53 * t76) - g(3) * (t52 * t51 + t53 * t74), t51 * t79 - g(2) * (t51 * t76 + t54 * t53) - g(3) * (-t51 * t74 + t52 * t53), t45, -g(1) * (t59 * t50 + t61 * t63 + t62) - g(2) * (pkin(4) * t75 + t69) - g(3) * (t50 * t74 - t54 * t72 + t70) + (-g(2) * (-t50 * t61 - pkin(2) + t72) - g(3) * pkin(4) * t58) * t52;];
U_reg = t1;
