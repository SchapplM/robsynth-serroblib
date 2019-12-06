% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:59
% EndTime: 2019-12-05 15:17:02
% DurationCPUTime: 0.73s
% Computational Cost: add. (567->106), mult. (1696->212), div. (0->0), fcn. (1607->8), ass. (0->70)
t42 = sin(qJ(3));
t43 = cos(qJ(4));
t79 = cos(qJ(5));
t56 = t79 * t43;
t40 = sin(qJ(5));
t41 = sin(qJ(4));
t76 = t40 * t41;
t48 = t56 - t76;
t19 = t48 * t42;
t81 = qJD(4) + qJD(5);
t80 = -pkin(7) - pkin(6);
t26 = t40 * t43 + t79 * t41;
t12 = t81 * t26;
t78 = t12 * t42;
t38 = sin(pkin(9));
t44 = cos(qJ(3));
t77 = t38 * t44;
t35 = t41 ^ 2;
t36 = t43 ^ 2;
t75 = t35 + t36;
t74 = qJD(3) * t38;
t73 = qJD(5) * t40;
t72 = t41 * qJD(4);
t71 = t42 * qJD(3);
t70 = t43 * qJD(4);
t69 = t44 * qJD(3);
t68 = -0.2e1 * pkin(3) * qJD(4);
t67 = pkin(4) * t72;
t66 = pkin(4) * t73;
t65 = t40 * t80;
t64 = t26 * t69;
t63 = t41 * t70;
t62 = t42 * t72;
t61 = t38 * t71;
t60 = t42 * t69;
t59 = t44 * t72;
t58 = t38 * t69;
t57 = t41 * t69;
t55 = t75 * t44;
t54 = t79 * qJD(5);
t53 = t41 * t65;
t52 = pkin(4) * t54;
t51 = t80 * t79;
t50 = t41 * t51;
t39 = cos(pkin(9));
t21 = -t39 * t41 + t43 * t77;
t49 = t39 * t43 + t41 * t77;
t8 = t79 * t21 - t40 * t49;
t47 = -t43 * t69 + t62;
t46 = t42 * t70 + t57;
t29 = t80 * t43;
t16 = -t79 * t29 + t53;
t13 = -qJD(4) * t21 + t41 * t61;
t14 = t49 * qJD(4) + t43 * t61;
t45 = -t13 * t41 - t14 * t43 + (-t21 * t41 + t43 * t49) * qJD(4);
t37 = t44 ^ 2;
t34 = -t43 * pkin(4) - pkin(3);
t30 = t42 ^ 2 * t74;
t27 = t38 ^ 2 * t60;
t18 = t26 * t42;
t15 = t40 * t29 + t50;
t11 = -qJD(4) * t56 - t43 * t54 + t81 * t76;
t7 = -t40 * t21 - t49 * t79;
t6 = -t16 * qJD(5) + (t43 * t51 - t53) * qJD(4);
t5 = -t29 * t73 - t81 * t50 - t65 * t70;
t4 = -t81 * t19 - t64;
t3 = t40 * t57 - t56 * t69 + t78;
t2 = -t8 * qJD(5) + t79 * t13 + t40 * t14;
t1 = -t40 * t13 + t79 * t14 + t21 * t73 + t49 * t54;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t13 * t49 - 0.2e1 * t21 * t14 + 0.2e1 * t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t8 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 + (-t37 * t38 + (t21 * t43 + t41 * t49) * t44) * qJD(3) + t45 * t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t19 - t2 * t18 - t8 * t3 - t37 * t74 + t7 * t4 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t75) * t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t18 * t4 - 0.2e1 * t19 * t3 - 0.2e1 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t61, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t38, t46 * t38, t45, -pkin(3) * t58 + t45 * pkin(6), 0, 0, 0, 0, 0, 0, (-t48 * t69 + t78) * t38, (-t11 * t42 + t64) * t38, -t1 * t48 + t7 * t11 - t8 * t12 - t2 * t26, -t1 * t16 + t2 * t15 - t8 * t5 + t7 * t6 + (pkin(4) * t62 + t34 * t69) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t69, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t71 - t59, t41 * t71 - t44 * t70, qJD(3) * t55, (-pkin(3) * t42 + pkin(6) * t55) * qJD(3), 0, 0, 0, 0, 0, 0, -t44 * t12 - t48 * t71, t44 * t11 + t26 * t71, -t18 * t11 - t19 * t12 - t4 * t26 - t3 * t48, -pkin(4) * t59 + t4 * t15 - t3 * t16 - t18 * t6 - t19 * t5 + t34 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t63, 0.2e1 * (-t35 + t36) * qJD(4), 0, -0.2e1 * t63, 0, 0, t41 * t68, t43 * t68, 0, 0, -0.2e1 * t26 * t11, -0.2e1 * t11 * t48 - 0.2e1 * t26 * t12, 0, -0.2e1 * t48 * t12, 0, 0, 0.2e1 * t34 * t12 - 0.2e1 * t48 * t67, -0.2e1 * t34 * t11 + 0.2e1 * t26 * t67, 0.2e1 * t15 * t11 - 0.2e1 * t16 * t12 - 0.2e1 * t6 * t26 - 0.2e1 * t48 * t5, 0.2e1 * t15 * t6 - 0.2e1 * t16 * t5 + 0.2e1 * t34 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, (t79 * t2 - t1 * t40 + (-t40 * t7 + t79 * t8) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t47, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, (t79 * t4 - t3 * t40 + (t18 * t40 + t79 * t19) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, -t72, 0, -pkin(6) * t70, pkin(6) * t72, 0, 0, 0, 0, -t11, 0, -t12, 0, t6, t5, (t79 * t11 - t12 * t40 + (t26 * t40 + t48 * t79) * qJD(5)) * pkin(4), (t79 * t6 - t40 * t5 + (-t15 * t40 + t79 * t16) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t66, -0.2e1 * t52, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t12, 0, t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t52, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
