% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:45
% EndTime: 2020-01-03 12:11:48
% DurationCPUTime: 0.72s
% Computational Cost: add. (526->83), mult. (1001->134), div. (0->0), fcn. (1058->6), ass. (0->66)
t55 = sin(qJ(2));
t75 = t55 * pkin(1);
t41 = pkin(7) + t75;
t54 = sin(qJ(3));
t51 = t54 ^ 2;
t57 = cos(qJ(3));
t52 = t57 ^ 2;
t65 = t51 + t52;
t67 = t65 * t41;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t29 = t53 * t54 - t56 * t57;
t84 = 0.2e1 * t29;
t31 = t53 * t57 + t54 * t56;
t83 = 0.2e1 * t31;
t44 = -pkin(3) * t57 - pkin(2);
t58 = cos(qJ(2));
t73 = t58 * pkin(1);
t33 = t44 - t73;
t82 = 0.2e1 * t33;
t81 = 0.2e1 * t44;
t80 = 0.2e1 * t57;
t25 = (-pkin(8) - t41) * t54;
t50 = t57 * pkin(8);
t71 = t57 * t41;
t26 = t50 + t71;
t14 = t56 * t25 - t26 * t53;
t63 = t31 * qJ(5);
t5 = t14 - t63;
t15 = t25 * t53 + t26 * t56;
t64 = t29 * qJ(5);
t6 = t15 - t64;
t79 = -t6 * t29 - t5 * t31;
t34 = (-pkin(8) - pkin(7)) * t54;
t74 = t57 * pkin(7);
t35 = t50 + t74;
t19 = t34 * t53 + t35 * t56;
t10 = t19 - t64;
t18 = t56 * t34 - t35 * t53;
t9 = t18 - t63;
t78 = -t10 * t29 - t9 * t31;
t77 = -t14 * t31 - t15 * t29;
t76 = t53 * pkin(3);
t48 = t56 * pkin(3);
t43 = -pkin(2) - t73;
t72 = pkin(2) - t43;
t70 = -t18 * t31 - t19 * t29;
t21 = pkin(4) * t29 + t44;
t20 = t21 - t73;
t69 = t20 + t21;
t68 = t33 + t44;
t66 = t65 * pkin(7);
t61 = pkin(3) ^ 2;
t59 = 0.2e1 * pkin(4);
t46 = t53 ^ 2 * t61;
t45 = -0.2e1 * t76;
t42 = t48 + pkin(4);
t38 = t54 * t80;
t28 = t31 ^ 2;
t27 = t29 ^ 2;
t24 = t31 * pkin(4);
t23 = t29 * t76;
t17 = -0.2e1 * t31 * t29;
t16 = -t31 * t48 - t23;
t13 = -t31 * t42 - t23;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t73, -0.2e1 * t75, 0, (t55 ^ 2 + t58 ^ 2) * pkin(1) ^ 2, t51, t38, 0, t52, 0, 0, -0.2e1 * t43 * t57, 0.2e1 * t43 * t54, 0.2e1 * t67, t41 ^ 2 * t65 + t43 ^ 2, t28, t17, 0, t27, 0, 0, t29 * t82, t31 * t82, 0.2e1 * t77, t14 ^ 2 + t15 ^ 2 + t33 ^ 2, t28, t17, 0, t27, 0, 0, t20 * t84, t20 * t83, 0.2e1 * t79, t20 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t73, -t75, 0, 0, t51, t38, 0, t52, 0, 0, t72 * t57, -t72 * t54, t66 + t67, -t43 * pkin(2) + pkin(7) * t67, t28, t17, 0, t27, 0, 0, t68 * t29, t68 * t31, t70 + t77, t14 * t18 + t15 * t19 + t33 * t44, t28, t17, 0, t27, 0, 0, t69 * t29, t69 * t31, t78 + t79, t10 * t6 + t20 * t21 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t51, t38, 0, t52, 0, 0, pkin(2) * t80, -0.2e1 * pkin(2) * t54, 0.2e1 * t66, pkin(7) ^ 2 * t65 + pkin(2) ^ 2, t28, t17, 0, t27, 0, 0, t29 * t81, t31 * t81, 0.2e1 * t70, t18 ^ 2 + t19 ^ 2 + t44 ^ 2, t28, t17, 0, t27, 0, 0, t21 * t84, t21 * t83, 0.2e1 * t78, t10 ^ 2 + t21 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t57, 0, -t54 * t41, -t71, 0, 0, 0, 0, t31, 0, -t29, 0, t14, -t15, t16, (t14 * t56 + t15 * t53) * pkin(3), 0, 0, t31, 0, -t29, 0, t5, -t6, t13, t42 * t5 + t6 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t57, 0, -t54 * pkin(7), -t74, 0, 0, 0, 0, t31, 0, -t29, 0, t18, -t19, t16, (t18 * t56 + t19 * t53) * pkin(3), 0, 0, t31, 0, -t29, 0, t9, -t10, t13, t10 * t76 + t42 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t48, t45, 0, t56 ^ 2 * t61 + t46, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, t45, 0, t42 ^ 2 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t29, 0, t14, -t15, 0, 0, 0, 0, t31, 0, -t29, 0, t5, -t6, -t24, t5 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t29, 0, t18, -t19, 0, 0, 0, 0, t31, 0, -t29, 0, t9, -t10, -t24, t9 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t48, -t76, 0, 0, 0, 0, 0, 0, 0, 1, t59 + t48, -t76, 0, t42 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t59, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t31, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t31, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
