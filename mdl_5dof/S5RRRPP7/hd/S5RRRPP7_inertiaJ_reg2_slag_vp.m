% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:54
% EndTime: 2019-12-31 21:05:57
% DurationCPUTime: 0.83s
% Computational Cost: add. (356->105), mult. (746->176), div. (0->0), fcn. (677->4), ass. (0->71)
t40 = sin(qJ(3));
t36 = t40 ^ 2;
t42 = cos(qJ(3));
t38 = t42 ^ 2;
t80 = t36 + t38;
t44 = pkin(3) + pkin(4);
t61 = t42 * qJ(4);
t79 = t44 * t40 - t61;
t62 = t40 * qJ(4);
t78 = t42 * t44 + t62;
t9 = pkin(2) + t78;
t77 = 0.2e1 * t9;
t52 = -t42 * pkin(3) - t62;
t14 = -pkin(2) + t52;
t76 = -0.2e1 * t14;
t41 = sin(qJ(2));
t75 = 0.2e1 * t41;
t74 = pkin(2) * t40;
t73 = pkin(2) * t42;
t37 = t41 ^ 2;
t72 = t37 * pkin(6);
t71 = t40 * pkin(7);
t70 = t41 * pkin(6);
t69 = t40 * t41;
t68 = t40 * t42;
t43 = cos(qJ(2));
t67 = t40 * t43;
t66 = t41 * t43;
t29 = t42 * t41;
t65 = t42 * t43;
t15 = -t43 * pkin(2) - t41 * pkin(7) - pkin(1);
t64 = pkin(6) * t67 - t42 * t15;
t7 = pkin(6) * t65 + t40 * t15;
t63 = t80 * pkin(7) ^ 2;
t60 = t42 * qJ(5);
t59 = t43 * qJ(4);
t58 = t40 * t66;
t57 = t37 * t68;
t56 = t41 * t65;
t35 = t43 * pkin(3);
t5 = t35 + t64;
t55 = -0.2e1 * t59 + t7;
t4 = -t59 + t7;
t54 = t4 * t42 + t5 * t40;
t53 = t40 * t64 + t7 * t42;
t51 = -pkin(3) * t40 + t61;
t50 = t41 * t60 - t5;
t49 = pkin(6) ^ 2;
t47 = qJ(4) ^ 2;
t46 = 0.2e1 * qJ(4);
t39 = t43 ^ 2;
t34 = t42 * pkin(7);
t32 = t37 * t49;
t28 = t38 * t37;
t27 = t36 * t37;
t26 = -0.2e1 * t68;
t24 = pkin(7) * t67;
t22 = qJ(5) * t69;
t21 = t40 * t29;
t20 = -0.2e1 * t56;
t19 = 0.2e1 * t57;
t18 = 0.2e1 * t58;
t17 = t34 - t60;
t16 = (pkin(7) - qJ(5)) * t40;
t13 = 0.2e1 * t80 * pkin(7);
t12 = (t36 - t38) * t41;
t8 = (pkin(6) - t51) * t41;
t3 = (-pkin(6) - t79) * t41;
t2 = t22 + t4;
t1 = t43 * pkin(4) - t50;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, 0.2e1 * t66, 0, t39, 0, 0, 0.2e1 * pkin(1) * t43, -0.2e1 * pkin(1) * t41, 0.2e1 * (t37 + t39) * pkin(6), pkin(1) ^ 2 + t39 * t49 + t32, t28, -0.2e1 * t57, t20, t27, t18, t39, 0.2e1 * t40 * t72 + 0.2e1 * t43 * t64, 0.2e1 * t42 * t72 + 0.2e1 * t7 * t43, (-t40 * t7 + t42 * t64) * t75, t64 ^ 2 + t7 ^ 2 + t32, t28, t20, t19, t39, -0.2e1 * t58, t27, 0.2e1 * t5 * t43 + 0.2e1 * t8 * t69, (-t4 * t40 + t42 * t5) * t75, -0.2e1 * t8 * t29 - 0.2e1 * t4 * t43, t4 ^ 2 + t5 ^ 2 + t8 ^ 2, t28, t19, 0.2e1 * t56, t27, t18, t39, 0.2e1 * t1 * t43 - 0.2e1 * t3 * t69, -0.2e1 * t2 * t43 + 0.2e1 * t3 * t29, (-t1 * t42 + t2 * t40) * t75, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t43, 0, -t70, -t43 * pkin(6), 0, 0, t21, -t12, -t67, -t21, -t65, 0, t24 + (-pkin(6) * t42 - t74) * t41, pkin(7) * t65 + (pkin(6) * t40 - t73) * t41, t53, -pkin(2) * t70 + t53 * pkin(7), t21, -t67, t12, 0, t65, -t21, t14 * t69 - t8 * t42 + t24, t54, -t8 * t40 + (-pkin(7) * t43 - t14 * t41) * t42, t54 * pkin(7) + t8 * t14, t21, t12, t67, -t21, -t65, 0, t16 * t43 + t3 * t42 - t9 * t69, -t17 * t43 + t9 * t29 + t3 * t40, (-t16 * t41 - t2) * t42 + (t17 * t41 - t1) * t40, t1 * t16 + t2 * t17 + t3 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, 0.2e1 * t68, 0, t38, 0, 0, 0.2e1 * t73, -0.2e1 * t74, t13, pkin(2) ^ 2 + t63, t36, 0, t26, 0, 0, t38, t42 * t76, t13, t40 * t76, t14 ^ 2 + t63, t36, t26, 0, t38, 0, 0, t42 * t77, t40 * t77, -0.2e1 * t16 * t40 - 0.2e1 * t17 * t42, t16 ^ 2 + t17 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t69, -t43, -t64, -t7, 0, 0, 0, t29, 0, -t43, t69, 0, -0.2e1 * t35 - t64, t52 * t41, t55, -t5 * pkin(3) + t4 * qJ(4), 0, 0, -t29, 0, -t69, -t43, (-pkin(4) - t44) * t43 + t50, t22 + t55, t78 * t41, t2 * qJ(4) - t1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t42, 0, -t71, -t34, 0, 0, 0, t40, 0, 0, -t42, 0, -t71, t51, t34, t51 * pkin(7), 0, 0, -t40, 0, t42, 0, -t16, t17, t79, t17 * qJ(4) - t16 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t46, pkin(3) ^ 2 + t47, 0, 0, 0, 0, 0, 1, 0.2e1 * t44, t46, 0, t44 ^ 2 + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t29, 0, t5, 0, 0, 0, 0, 0, 0, t43, 0, -t29, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t71, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -1, 0, 0, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t29, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t40, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
