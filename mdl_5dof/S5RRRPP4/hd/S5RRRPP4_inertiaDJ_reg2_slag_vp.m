% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:46
% EndTime: 2019-12-31 20:55:49
% DurationCPUTime: 0.89s
% Computational Cost: add. (1566->118), mult. (3704->212), div. (0->0), fcn. (3261->6), ass. (0->69)
t56 = sin(qJ(3));
t58 = cos(qJ(3));
t57 = sin(qJ(2));
t87 = -pkin(7) - pkin(6);
t73 = t87 * t57;
t59 = cos(qJ(2));
t89 = t87 * t59;
t90 = -t56 * t73 + t58 * t89;
t25 = t56 * t89 + t58 * t73;
t88 = qJD(2) + qJD(3);
t60 = 2 * qJD(5);
t84 = t56 * t57;
t40 = -t58 * t59 + t84;
t18 = -t40 * qJ(4) - t90;
t55 = sin(pkin(8));
t41 = t56 * t59 + t58 * t57;
t61 = -t41 * qJ(4) + t25;
t82 = cos(pkin(8));
t11 = t55 * t18 - t82 * t61;
t69 = t82 * t56;
t83 = pkin(2) * qJD(3);
t32 = (t55 * t58 + t69) * t83;
t86 = t11 * t32;
t22 = -t55 * t40 + t82 * t41;
t85 = t32 * t22;
t52 = t58 * pkin(2) + pkin(3);
t35 = pkin(2) * t69 + t55 * t52;
t81 = qJD(3) * t58;
t80 = t57 * qJD(2);
t79 = t59 * qJD(2);
t23 = -t58 * t79 - t59 * t81 + t88 * t84;
t24 = t88 * t41;
t13 = -t55 * t23 + t82 * t24;
t21 = t82 * t40 + t55 * t41;
t78 = 0.2e1 * t21 * t13;
t77 = -0.2e1 * pkin(1) * qJD(2);
t54 = pkin(2) * t80;
t76 = t56 * t83;
t75 = pkin(2) * t81;
t71 = t57 * t79;
t53 = -t59 * pkin(2) - pkin(1);
t12 = t82 * t18 + t55 * t61;
t62 = t90 * qJD(2);
t16 = qJD(3) * t90 + t62;
t64 = t23 * qJ(4) - t41 * qJD(4);
t66 = qJD(3) * t89;
t67 = qJD(3) * t73;
t15 = -t25 * qJD(2) - t56 * t66 - t58 * t67;
t8 = -t24 * qJ(4) - t40 * qJD(4) - t15;
t3 = t55 * t8 - t82 * (t16 + t64);
t4 = t82 * t8 + (-t56 * t67 + t58 * t66 + t62 + t64) * t55;
t70 = t11 * t3 + t12 * t4;
t19 = t24 * pkin(3) + t54;
t14 = -t82 * t23 - t55 * t24;
t65 = t22 * t13 + t14 * t21;
t28 = t40 * pkin(3) + t53;
t33 = -t55 * t76 + t82 * t75;
t34 = -t55 * t56 * pkin(2) + t82 * t52;
t63 = 0.2e1 * t11 * t14 - 0.2e1 * t12 * t13 - 0.2e1 * t4 * t21 + 0.2e1 * t3 * t22;
t50 = -t82 * pkin(3) - pkin(4);
t49 = t55 * pkin(3) + qJ(5);
t31 = -pkin(4) - t34;
t30 = qJ(5) + t35;
t29 = qJD(5) + t33;
t27 = 0.2e1 * t32;
t10 = t21 * pkin(4) - t22 * qJ(5) + t28;
t9 = 0.2e1 * t22 * t14;
t5 = t13 * pkin(4) - t14 * qJ(5) - t22 * qJD(5) + t19;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t71, 0.2e1 * (-t57 ^ 2 + t59 ^ 2) * qJD(2), 0, -0.2e1 * t71, 0, 0, t57 * t77, t59 * t77, 0, 0, -0.2e1 * t41 * t23, 0.2e1 * t23 * t40 - 0.2e1 * t41 * t24, 0, 0.2e1 * t40 * t24, 0, 0, 0.2e1 * t53 * t24 + 0.2e1 * t40 * t54, -0.2e1 * t53 * t23 + 0.2e1 * t41 * t54, 0.2e1 * t15 * t40 - 0.2e1 * t16 * t41 + 0.2e1 * t25 * t23 + 0.2e1 * t24 * t90, 0.2e1 * t15 * t90 + 0.2e1 * t25 * t16 + 0.2e1 * t53 * t54, t9, -0.2e1 * t65, 0, t78, 0, 0, 0.2e1 * t28 * t13 + 0.2e1 * t19 * t21, 0.2e1 * t28 * t14 + 0.2e1 * t19 * t22, t63, 0.2e1 * t28 * t19 + 0.2e1 * t70, t9, 0, 0.2e1 * t65, 0, 0, t78, 0.2e1 * t10 * t13 + 0.2e1 * t5 * t21, t63, -0.2e1 * t10 * t14 - 0.2e1 * t5 * t22, 0.2e1 * t10 * t5 + 0.2e1 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, -t80, 0, -pkin(6) * t79, pkin(6) * t80, 0, 0, 0, 0, -t23, 0, -t24, 0, t16, t15, (t23 * t58 - t24 * t56 + (-t40 * t58 + t41 * t56) * qJD(3)) * pkin(2), (-t15 * t56 + t16 * t58 + (-t25 * t56 - t58 * t90) * qJD(3)) * pkin(2), 0, 0, t14, 0, -t13, 0, -t3, -t4, -t35 * t13 - t34 * t14 - t33 * t21 + t85, t12 * t33 - t3 * t34 + t4 * t35 + t86, 0, t14, 0, 0, t13, 0, -t3, -t30 * t13 + t31 * t14 - t29 * t21 + t85, t4, t12 * t29 + t3 * t31 + t4 * t30 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t76, -0.2e1 * t75, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -0.2e1 * t33, 0, -0.2e1 * t34 * t32 + 0.2e1 * t35 * t33, 0, 0, 0, 0, 0, 0, -t27, 0, 0.2e1 * t29, 0.2e1 * t30 * t29 + 0.2e1 * t31 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t24, 0, t16, t15, 0, 0, 0, 0, t14, 0, -t13, 0, -t3, -t4, (-t13 * t55 - t82 * t14) * pkin(3), (-t82 * t3 + t4 * t55) * pkin(3), 0, t14, 0, 0, t13, 0, -t3, -qJD(5) * t21 - t49 * t13 + t50 * t14, t4, t12 * qJD(5) + t3 * t50 + t4 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, (-t82 * t32 + t33 * t55) * pkin(3), 0, 0, 0, 0, 0, 0, -t32, 0, t60 + t33, t30 * qJD(5) + t29 * t49 + t32 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t49 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, t19, 0, 0, 0, 0, 0, 0, t13, 0, -t14, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
