% Calculate inertial parameters regressor of potential energy for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPPRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:36
% EndTime: 2019-12-05 14:59:37
% DurationCPUTime: 0.17s
% Computational Cost: add. (138->67), mult. (302->98), div. (0->0), fcn. (348->10), ass. (0->43)
t84 = g(3) * qJ(1);
t55 = sin(pkin(9));
t56 = sin(pkin(8));
t83 = t55 * t56;
t58 = cos(pkin(9));
t82 = t56 * t58;
t62 = sin(qJ(4));
t81 = t56 * t62;
t64 = cos(qJ(4));
t80 = t56 * t64;
t57 = sin(pkin(7));
t59 = cos(pkin(8));
t79 = t57 * t59;
t60 = cos(pkin(7));
t78 = t60 * t55;
t77 = t60 * t58;
t76 = t60 * pkin(1) + t57 * qJ(2);
t75 = qJ(3) * t56;
t74 = t57 * pkin(1) - t60 * qJ(2);
t73 = t76 + (pkin(2) * t59 + t75) * t60;
t72 = g(1) * t60 + g(2) * t57;
t71 = t56 * pkin(2) - t59 * qJ(3) + qJ(1);
t70 = pkin(2) * t79 + t57 * t75 + t74;
t69 = pkin(3) * t82 + pkin(5) * t83 + t71;
t36 = t58 * t79 - t78;
t28 = t36 * t62 - t57 * t80;
t38 = t57 * t55 + t59 * t77;
t30 = t38 * t62 - t60 * t80;
t39 = t58 * t81 + t59 * t64;
t68 = g(1) * t30 + g(2) * t28 + g(3) * t39;
t37 = -t57 * t58 + t59 * t78;
t67 = t38 * pkin(3) + t37 * pkin(5) + t73;
t35 = t55 * t79 + t77;
t66 = g(1) * t37 + g(2) * t35 + g(3) * t83;
t65 = t36 * pkin(3) + t35 * pkin(5) + t70;
t63 = cos(qJ(5));
t61 = sin(qJ(5));
t41 = g(1) * t57 - g(2) * t60;
t40 = t58 * t80 - t59 * t62;
t32 = -g(3) * t59 + t72 * t56;
t31 = t38 * t64 + t60 * t81;
t29 = t36 * t64 + t57 * t81;
t1 = [0, 0, 0, 0, 0, 0, -t72, t41, -g(3), -t84, 0, 0, 0, 0, 0, 0, -g(3) * t56 - t72 * t59, t32, -t41, -g(1) * t76 - g(2) * t74 - t84, 0, 0, 0, 0, 0, 0, -g(1) * t38 - g(2) * t36 - g(3) * t82, t66, -t32, -g(1) * t73 - g(2) * t70 - g(3) * t71, 0, 0, 0, 0, 0, 0, -g(1) * t31 - g(2) * t29 - g(3) * t40, t68, -t66, -g(1) * t67 - g(2) * t65 - g(3) * t69, 0, 0, 0, 0, 0, 0, -g(1) * (t31 * t63 + t37 * t61) - g(2) * (t29 * t63 + t35 * t61) - g(3) * (t40 * t63 + t61 * t83), -g(1) * (-t31 * t61 + t37 * t63) - g(2) * (-t29 * t61 + t35 * t63) - g(3) * (-t40 * t61 + t63 * t83), -t68, -g(1) * (t31 * pkin(4) + t30 * pkin(6) + t67) - g(2) * (t29 * pkin(4) + t28 * pkin(6) + t65) - g(3) * (t40 * pkin(4) + t39 * pkin(6) + t69);];
U_reg = t1;
