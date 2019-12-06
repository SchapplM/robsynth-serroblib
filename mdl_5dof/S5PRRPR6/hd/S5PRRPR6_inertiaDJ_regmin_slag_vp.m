% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:47
% EndTime: 2019-12-05 16:32:51
% DurationCPUTime: 0.82s
% Computational Cost: add. (541->130), mult. (1647->271), div. (0->0), fcn. (1529->10), ass. (0->82)
t49 = sin(pkin(10));
t51 = cos(pkin(10));
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t91 = -t53 * t49 + t56 * t51;
t45 = -t51 * pkin(4) - pkin(3);
t90 = 0.2e1 * t45;
t54 = sin(qJ(3));
t89 = t49 * t54;
t57 = cos(qJ(3));
t88 = t49 * t57;
t50 = sin(pkin(5));
t55 = sin(qJ(2));
t87 = t50 * t55;
t58 = cos(qJ(2));
t86 = t50 * t58;
t85 = t51 * t54;
t84 = t51 * t57;
t81 = pkin(8) + qJ(4);
t27 = -t54 * qJD(4) + (pkin(3) * t54 - qJ(4) * t57) * qJD(3);
t77 = t54 * qJD(3);
t73 = pkin(7) * t77;
t14 = t51 * t27 + t49 * t73;
t66 = -t57 * pkin(3) - t54 * qJ(4);
t38 = -pkin(2) + t66;
t43 = pkin(7) * t84;
t22 = t49 * t38 + t43;
t80 = qJD(2) * t55;
t79 = qJD(4) * t57;
t78 = qJD(5) * t54;
t76 = t57 * qJD(3);
t75 = pkin(7) * t88;
t74 = -0.2e1 * pkin(2) * qJD(3);
t46 = pkin(7) * t76;
t72 = t49 * t76;
t71 = t50 * t80;
t70 = qJD(2) * t86;
t69 = t54 * t76;
t68 = 0.2e1 * (t49 ^ 2 + t51 ^ 2) * qJD(4);
t52 = cos(pkin(5));
t30 = -t52 * t57 + t54 * t87;
t18 = -qJD(3) * t30 + t57 * t70;
t8 = -t18 * t49 + t51 * t71;
t9 = t18 * t51 + t49 * t71;
t67 = -t8 * t49 + t9 * t51;
t34 = t51 * t38;
t13 = -pkin(8) * t85 + t34 + (-pkin(7) * t49 - pkin(4)) * t57;
t20 = -pkin(8) * t89 + t22;
t65 = t56 * t13 - t53 * t20;
t64 = t53 * t13 + t56 * t20;
t23 = t49 * t27;
t15 = -t51 * t73 + t23;
t63 = -t14 * t49 + t15 * t51;
t31 = t52 * t54 + t57 * t87;
t16 = -t31 * t49 - t51 * t86;
t17 = t31 * t51 - t49 * t86;
t62 = t56 * t16 - t53 * t17;
t61 = t53 * t16 + t56 * t17;
t40 = t81 * t49;
t41 = t81 * t51;
t60 = -t56 * t40 - t53 * t41;
t59 = -t53 * t40 + t56 * t41;
t36 = t56 * t49 + t53 * t51;
t37 = (pkin(4) * t49 + pkin(7)) * t54;
t32 = pkin(4) * t72 + t46;
t29 = t36 * qJD(5);
t28 = t91 * qJD(5);
t26 = t91 * t54;
t25 = t36 * t54;
t21 = t34 - t75;
t19 = qJD(3) * t31 + t54 * t70;
t12 = t23 + (-pkin(7) * t85 - pkin(8) * t88) * qJD(3);
t11 = t36 * t76 + t91 * t78;
t10 = -t36 * t78 + t76 * t91;
t7 = (pkin(4) * t54 - pkin(8) * t84) * qJD(3) + t14;
t6 = -qJD(4) * t36 - qJD(5) * t59;
t5 = -qJD(4) * t91 - qJD(5) * t60;
t4 = -qJD(5) * t64 - t53 * t12 + t56 * t7;
t3 = -qJD(5) * t65 - t56 * t12 - t53 * t7;
t2 = -qJD(5) * t61 - t53 * t9 + t56 * t8;
t1 = -qJD(5) * t62 - t53 * t8 - t56 * t9;
t24 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t16 * t8 + 0.2e1 * t17 * t9 + 0.2e1 * t30 * t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t71, -t70, 0, 0, 0, 0, 0, (-t57 * t80 - t58 * t77) * t50, (t54 * t80 - t58 * t76) * t50, t19 * t89 - t8 * t57 + (t16 * t54 + t30 * t88) * qJD(3), t19 * t85 + t9 * t57 + (-t17 * t54 + t30 * t84) * qJD(3), (-t49 * t9 - t51 * t8) * t54 + (-t16 * t51 - t17 * t49) * t76, t16 * t14 + t17 * t15 + t8 * t21 + t9 * t22 + (t19 * t54 + t30 * t76) * pkin(7), 0, 0, 0, 0, 0, t30 * t11 + t19 * t25 - t2 * t57 + t62 * t77, -t1 * t57 + t30 * t10 + t19 * t26 - t61 * t77; 0, 0, 0, 0, 0.2e1 * t69, 0.2e1 * (-t54 ^ 2 + t57 ^ 2) * qJD(3), 0, 0, 0, t54 * t74, t57 * t74, -0.2e1 * t14 * t57 + 0.2e1 * (t21 + 0.2e1 * t75) * t77, 0.2e1 * t15 * t57 + 0.2e1 * (-t22 + 0.2e1 * t43) * t77, 0.2e1 * (-t14 * t51 - t15 * t49) * t54 + 0.2e1 * (-t21 * t51 - t22 * t49) * t76, 0.2e1 * pkin(7) ^ 2 * t69 + 0.2e1 * t21 * t14 + 0.2e1 * t22 * t15, 0.2e1 * t26 * t10, -0.2e1 * t10 * t25 - 0.2e1 * t26 * t11, -0.2e1 * t10 * t57 + 0.2e1 * t26 * t77, 0.2e1 * t11 * t57 - 0.2e1 * t25 * t77, -0.2e1 * t69, 0.2e1 * t37 * t11 + 0.2e1 * t32 * t25 - 0.2e1 * t4 * t57 + 0.2e1 * t65 * t77, 0.2e1 * t37 * t10 + 0.2e1 * t32 * t26 - 0.2e1 * t3 * t57 - 0.2e1 * t64 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t18, -t19 * t51, t19 * t49, t67, -t19 * pkin(3) + (-t16 * t49 + t17 * t51) * qJD(4) + t67 * qJ(4), 0, 0, 0, 0, 0, -t19 * t91 + t30 * t29, t19 * t36 + t30 * t28; 0, 0, 0, 0, 0, 0, t76, -t77, 0, -t46, t73, t49 * t79 + (t49 * t66 - t43) * qJD(3), t51 * t79 + (t51 * t66 + t75) * qJD(3), t63, -pkin(3) * t46 + (-t21 * t49 + t22 * t51) * qJD(4) + t63 * qJ(4), t10 * t36 + t26 * t28, t10 * t91 - t36 * t11 - t28 * t25 - t26 * t29, -t28 * t57 + t36 * t77, t29 * t57 + t77 * t91, 0, t45 * t11 + t37 * t29 - t32 * t91 - t6 * t57 + t60 * t77, t45 * t10 + t37 * t28 + t32 * t36 - t5 * t57 - t59 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, qJ(4) * t68, 0.2e1 * t36 * t28, 0.2e1 * t28 * t91 - 0.2e1 * t36 * t29, 0, 0, 0, t29 * t90, t28 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t51 * t76, 0, t46, 0, 0, 0, 0, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t77, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t24;
