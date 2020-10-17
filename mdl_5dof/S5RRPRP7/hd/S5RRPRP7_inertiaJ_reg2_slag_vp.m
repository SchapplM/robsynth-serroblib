% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:31
% EndTime: 2019-12-31 20:01:34
% DurationCPUTime: 0.73s
% Computational Cost: add. (597->80), mult. (1171->156), div. (0->0), fcn. (1292->6), ass. (0->67)
t45 = sin(qJ(4));
t39 = t45 ^ 2;
t47 = cos(qJ(4));
t41 = t47 ^ 2;
t33 = t39 + t41;
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t28 = t43 * t48 + t44 * t46;
t79 = -0.2e1 * t28;
t66 = -qJ(3) - pkin(6);
t30 = t66 * t48;
t59 = t66 * t46;
t13 = -t43 * t30 - t44 * t59;
t78 = t13 ^ 2;
t26 = t43 * t46 - t44 * t48;
t24 = t26 ^ 2;
t38 = -t48 * pkin(2) - pkin(1);
t77 = 0.2e1 * t38;
t76 = -0.2e1 * t47;
t75 = 0.2e1 * t48;
t74 = t26 * pkin(4);
t73 = t43 * pkin(2);
t72 = t44 * pkin(2);
t15 = -t44 * t30 + t43 * t59;
t8 = t26 * pkin(3) - t28 * pkin(7) + t38;
t4 = t47 * t15 + t45 * t8;
t36 = pkin(7) + t73;
t71 = t26 * t36;
t19 = t45 * t26;
t70 = t45 * t28;
t69 = t45 * t36;
t68 = t45 * t47;
t21 = t47 * t26;
t22 = t47 * t28;
t67 = t47 * t36;
t65 = t33 * t36 ^ 2;
t40 = t46 ^ 2;
t42 = t48 ^ 2;
t64 = t40 + t42;
t63 = t26 * qJ(5);
t62 = t26 * t70;
t25 = t28 ^ 2;
t61 = t25 * t68;
t37 = -pkin(3) - t72;
t60 = t45 * t15 - t47 * t8;
t1 = t63 + t4;
t2 = t60 - t74;
t58 = t1 * t47 + t2 * t45;
t57 = t1 * t45 - t2 * t47;
t56 = t4 * t45 - t47 * t60;
t55 = t4 * t47 + t45 * t60;
t54 = t47 * pkin(4) + t45 * qJ(5);
t53 = pkin(4) * t45 - t47 * qJ(5);
t23 = t37 - t54;
t52 = t23 * t28 - t71;
t51 = t28 * t37 - t71;
t20 = t41 * t25;
t18 = t39 * t25;
t17 = t45 * t22;
t16 = 0.2e1 * t33 * t36;
t11 = 0.2e1 * t26 * t22;
t10 = t33 * t28;
t9 = (t39 - t41) * t28;
t5 = t53 * t28 + t13;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t40, t46 * t75, 0, t42, 0, 0, pkin(1) * t75, -0.2e1 * pkin(1) * t46, 0.2e1 * t64 * pkin(6), t64 * pkin(6) ^ 2 + pkin(1) ^ 2, t25, t26 * t79, 0, t24, 0, 0, t26 * t77, t28 * t77, 0.2e1 * t13 * t28 - 0.2e1 * t15 * t26, t15 ^ 2 + t38 ^ 2 + t78, t20, -0.2e1 * t61, t11, t18, -0.2e1 * t62, t24, 0.2e1 * t13 * t70 - 0.2e1 * t26 * t60, 0.2e1 * t13 * t22 - 0.2e1 * t4 * t26, t56 * t79, t4 ^ 2 + t60 ^ 2 + t78, t20, t11, 0.2e1 * t61, t24, 0.2e1 * t62, t18, -0.2e1 * t2 * t26 + 0.2e1 * t5 * t70, t57 * t79, 0.2e1 * t1 * t26 - 0.2e1 * t5 * t22, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t48, 0, -t46 * pkin(6), -t48 * pkin(6), 0, 0, 0, 0, t28, 0, -t26, 0, -t13, -t15, (-t26 * t43 - t28 * t44) * pkin(2), (-t13 * t44 + t15 * t43) * pkin(2), t17, -t9, t19, -t17, t21, 0, -t13 * t47 + t51 * t45, t13 * t45 + t51 * t47, t55, t13 * t37 + t55 * t36, t17, t19, t9, 0, -t21, -t17, t52 * t45 - t5 * t47, t58, -t5 * t45 - t52 * t47, t5 * t23 + t58 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t72, -0.2e1 * t73, 0, (t43 ^ 2 + t44 ^ 2) * pkin(2) ^ 2, t39, 0.2e1 * t68, 0, t41, 0, 0, t37 * t76, 0.2e1 * t37 * t45, t16, t37 ^ 2 + t65, t39, 0, -0.2e1 * t68, 0, 0, t41, t23 * t76, t16, -0.2e1 * t23 * t45, t23 ^ 2 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t28, 0, t38, 0, 0, 0, 0, 0, 0, t21, -t19, -t10, t56, 0, 0, 0, 0, 0, 0, t21, -t10, t19, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t70, t26, -t60, -t4, 0, 0, 0, t22, 0, t26, t70, 0, -t60 + 0.2e1 * t74, -t54 * t28, 0.2e1 * t63 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t47, 0, -t69, -t67, 0, 0, 0, t45, 0, 0, -t47, 0, -t69, -t53, t67, -t53 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t45, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t45, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t22, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
