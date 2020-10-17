% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP13_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:37
% EndTime: 2019-12-31 18:59:40
% DurationCPUTime: 0.71s
% Computational Cost: add. (272->74), mult. (535->132), div. (0->0), fcn. (478->4), ass. (0->65)
t37 = sin(qJ(4));
t32 = t37 ^ 2;
t39 = cos(qJ(4));
t34 = t39 ^ 2;
t53 = t32 + t34;
t71 = 2 * pkin(4);
t70 = -0.2e1 * t37;
t69 = 0.2e1 * t39;
t40 = cos(qJ(3));
t68 = 0.2e1 * t40;
t67 = 2 * qJ(2);
t38 = sin(qJ(3));
t66 = pkin(7) * t38;
t65 = t37 * pkin(7);
t64 = t39 * pkin(7);
t63 = t40 * pkin(3);
t15 = t38 * pkin(3) - t40 * pkin(7) + qJ(2);
t41 = -pkin(1) - pkin(6);
t59 = t39 * t41;
t4 = t37 * t15 + t38 * t59;
t62 = t37 * t39;
t25 = t37 * t40;
t61 = t37 * t41;
t60 = t38 * t41;
t28 = t39 * t40;
t45 = -t39 * pkin(4) - t37 * qJ(5);
t16 = -pkin(3) + t45;
t58 = t40 * t16;
t57 = t40 * t38;
t56 = t40 * t41;
t55 = t53 * t66;
t54 = t53 * pkin(7) ^ 2;
t33 = t38 ^ 2;
t35 = t40 ^ 2;
t20 = t33 + t35;
t52 = t38 * qJ(5);
t51 = t37 * t57;
t50 = t35 * t62;
t49 = -t63 - t66;
t1 = t52 + t4;
t8 = t39 * t15;
t2 = -t8 + (-pkin(4) + t61) * t38;
t48 = t1 * t39 + t2 * t37;
t3 = -t37 * t60 + t8;
t47 = -t3 * t37 + t4 * t39;
t46 = -t58 + t66;
t44 = -pkin(4) * t37 + t39 * qJ(5);
t42 = qJ(2) ^ 2;
t36 = t41 ^ 2;
t29 = t35 * t36;
t27 = t39 * t38;
t26 = t34 * t35;
t24 = t37 * t38;
t23 = t32 * t35;
t18 = t37 * t28;
t17 = t57 * t69;
t14 = 0.2e1 * t53 * pkin(7);
t13 = t20 * t41;
t12 = t20 * t39;
t11 = t53 * t38;
t10 = t20 * t37;
t9 = (t32 - t34) * t40;
t6 = t53 * t33 + t35;
t5 = (-t41 - t44) * t40;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t67, pkin(1) ^ 2 + t42, t35, -0.2e1 * t57, 0, t33, 0, 0, t38 * t67, t40 * t67, -0.2e1 * t13, t33 * t36 + t29 + t42, t26, -0.2e1 * t50, t17, t23, -0.2e1 * t51, t33, 0.2e1 * t3 * t38 - 0.2e1 * t35 * t61, -0.2e1 * t35 * t59 - 0.2e1 * t4 * t38, (-t3 * t39 - t37 * t4) * t68, t3 ^ 2 + t4 ^ 2 + t29, t26, t17, 0.2e1 * t50, t33, 0.2e1 * t51, t23, -0.2e1 * t2 * t38 + 0.2e1 * t5 * t25, (-t1 * t37 + t2 * t39) * t68, 0.2e1 * t1 * t38 - 0.2e1 * t5 * t28, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t20, t13, 0, 0, 0, 0, 0, 0, -t10, -t12, 0, t35 * t41 + t47 * t38, 0, 0, 0, 0, 0, 0, -t10, 0, t12, t48 * t38 - t5 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t38, 0, t56, -t60, 0, 0, t18, -t9, t24, -t18, t27, 0, t49 * t37 + t39 * t56, -t37 * t56 + t49 * t39, t47, pkin(3) * t56 + t47 * pkin(7), t18, t24, t9, 0, -t27, -t18, -t46 * t37 - t5 * t39, t48, -t5 * t37 + t46 * t39, t48 * pkin(7) + t5 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t38, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t25, t11, t55 + t63, 0, 0, 0, 0, 0, 0, t28, t11, t25, t55 - t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t32, 0.2e1 * t62, 0, t34, 0, 0, pkin(3) * t69, pkin(3) * t70, t14, pkin(3) ^ 2 + t54, t32, 0, -0.2e1 * t62, 0, 0, t34, -0.2e1 * t16 * t39, t14, t16 * t70, t16 ^ 2 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t25, t38, t3, -t4, 0, 0, 0, t28, 0, t38, t25, 0, t8 + (t71 - t61) * t38, t45 * t40, 0.2e1 * t52 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t27, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, t27, t44 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t39, 0, -t65, -t64, 0, 0, 0, t37, 0, 0, -t39, 0, -t65, t44, t64, t44 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t71, 0, 0.2e1 * qJ(5), (pkin(4) ^ 2) + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t28, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
