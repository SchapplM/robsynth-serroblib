% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:40
% EndTime: 2019-12-31 19:01:42
% DurationCPUTime: 0.51s
% Computational Cost: add. (644->96), mult. (1517->175), div. (0->0), fcn. (1353->8), ass. (0->74)
t52 = cos(qJ(5));
t47 = t52 ^ 2;
t49 = sin(qJ(5));
t79 = t49 ^ 2 - t47;
t64 = t79 * qJD(5);
t90 = qJD(3) + qJD(4);
t51 = sin(qJ(3));
t40 = sin(pkin(9)) * pkin(1) + pkin(6);
t89 = pkin(7) + t40;
t27 = t89 * t51;
t53 = cos(qJ(3));
t28 = t89 * t53;
t50 = sin(qJ(4));
t87 = cos(qJ(4));
t19 = -t50 * t27 + t87 * t28;
t66 = t87 * qJD(3);
t56 = t89 * t66;
t62 = qJD(3) * t50 * t89;
t6 = t19 * qJD(4) - t51 * t62 + t53 * t56;
t4 = t6 * t49;
t18 = t87 * t27 + t50 * t28;
t45 = qJD(5) * t52;
t88 = t18 * t45 + t4;
t65 = t87 * qJD(4);
t84 = t50 * t51;
t21 = t90 * t84 + (-t66 - t65) * t53;
t32 = t50 * t53 + t87 * t51;
t86 = t32 * t21;
t22 = t90 * t32;
t85 = t49 * t22;
t83 = t52 * t21;
t82 = t52 * t22;
t31 = -t87 * t53 + t84;
t81 = -t31 * t83 + t32 * t82;
t44 = -t87 * pkin(3) - pkin(4);
t78 = qJD(4) * t50;
t69 = pkin(3) * t78;
t80 = t44 * t45 + t49 * t69;
t77 = qJD(5) * t49;
t76 = t51 * qJD(3);
t75 = t53 * qJD(3);
t74 = t49 * t83;
t73 = 0.2e1 * t75;
t72 = pkin(4) * t77;
t71 = pkin(4) * t45;
t70 = pkin(3) * t76;
t68 = t32 * t77;
t67 = t49 * t45;
t41 = -cos(pkin(9)) * pkin(1) - pkin(2);
t63 = pkin(3) * t65;
t35 = -t53 * pkin(3) + t41;
t17 = t31 * pkin(4) - t32 * pkin(8) + t35;
t61 = t52 * t17 - t49 * t19;
t60 = t49 * t17 + t52 * t19;
t59 = t21 * t31 - t32 * t22;
t43 = t50 * pkin(3) + pkin(8);
t58 = t31 * t43 - t32 * t44;
t57 = t44 * t77 - t52 * t69;
t55 = -t49 * t21 + t32 * t45;
t10 = t68 + t83;
t9 = t31 * t77 - t82;
t54 = -t21 * t44 - t22 * t43 + (-t87 * t31 + t32 * t50) * qJD(4) * pkin(3);
t37 = 0.2e1 * t67;
t30 = -0.2e1 * t64;
t29 = t32 ^ 2;
t15 = t18 * t77;
t11 = t31 * t45 + t85;
t8 = t22 * pkin(4) + t21 * pkin(8) + t70;
t7 = -t32 * t64 - t74;
t5 = t27 * t65 + t28 * t78 + t51 * t56 + t53 * t62;
t3 = t79 * t21 - 0.4e1 * t32 * t67;
t2 = -t60 * qJD(5) + t49 * t5 + t52 * t8;
t1 = -t61 * qJD(5) - t49 * t8 + t52 * t5;
t12 = [0, 0, 0, 0, t51 * t73, 0.2e1 * (-t51 ^ 2 + t53 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t41 * t76, t41 * t73, -0.2e1 * t86, 0.2e1 * t59, 0, 0, 0, 0.2e1 * t35 * t22 + 0.2e1 * t31 * t70, -0.2e1 * t35 * t21 + 0.2e1 * t32 * t70, -0.2e1 * t29 * t67 - 0.2e1 * t47 * t86, 0.2e1 * t29 * t64 + 0.4e1 * t32 * t74, -0.2e1 * t31 * t68 + 0.2e1 * t81, -0.2e1 * t31 * t55 - 0.2e1 * t32 * t85, 0.2e1 * t31 * t22, 0.2e1 * t18 * t55 + 0.2e1 * t2 * t31 + 0.2e1 * t22 * t61 + 0.2e1 * t32 * t4, 0.2e1 * t6 * t52 * t32 + 0.2e1 * t1 * t31 - 0.2e1 * t10 * t18 - 0.2e1 * t60 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t59 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t75, -t76, 0, -t40 * t75, t40 * t76, 0, 0, -t21, -t22, 0, -t6, t5, t7, t3, t11, -t9, 0, t15 + (-qJD(5) * t58 - t6) * t52 + t54 * t49, t52 * t54 + t58 * t77 + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75, 0, 0, 0, 0, 0, -t22, t21, 0, 0, 0, 0, 0, t9, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t69, -0.2e1 * t63, t37, t30, 0, 0, 0, 0.2e1 * t57, 0.2e1 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, 0, -t6, t5, t7, t3, t11, -t9, 0, t15 + (pkin(4) * t21 - pkin(8) * t22) * t49 + (-t6 + (-pkin(4) * t32 - pkin(8) * t31) * qJD(5)) * t52, pkin(4) * t10 + pkin(8) * t9 + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, 0, 0, 0, 0, 0, t9, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t63, t37, t30, 0, 0, 0, t57 - t72, -t71 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t30, 0, 0, 0, -0.2e1 * t72, -0.2e1 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t55, t22, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t77, 0, -t43 * t45 - t49 * t63, t43 * t77 - t52 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t77, 0, -pkin(8) * t45, pkin(8) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
