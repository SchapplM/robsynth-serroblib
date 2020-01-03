% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:11
% EndTime: 2019-12-31 18:22:13
% DurationCPUTime: 0.58s
% Computational Cost: add. (464->96), mult. (1206->198), div. (0->0), fcn. (1024->8), ass. (0->70)
t47 = sin(pkin(9));
t48 = cos(pkin(9));
t50 = sin(qJ(5));
t52 = cos(qJ(5));
t85 = -t50 * t47 + t52 * t48;
t26 = t85 * qJD(5);
t41 = -t48 * pkin(4) - pkin(3);
t84 = 0.2e1 * t41;
t51 = sin(qJ(3));
t83 = pkin(3) * t51;
t53 = cos(qJ(3));
t82 = t47 * t53;
t81 = t48 * t51;
t80 = t48 * t53;
t77 = pkin(7) + qJ(4);
t25 = -t51 * qJD(4) + (-qJ(4) * t53 + t83) * qJD(3);
t39 = sin(pkin(8)) * pkin(1) + pkin(6);
t42 = t51 * qJD(3);
t66 = t39 * t42;
t12 = t48 * t25 + t47 * t66;
t40 = -cos(pkin(8)) * pkin(1) - pkin(2);
t59 = -t53 * pkin(3) - t51 * qJ(4);
t28 = t40 + t59;
t32 = t39 * t80;
t15 = t47 * t28 + t32;
t76 = t47 ^ 2 + t48 ^ 2;
t75 = t51 ^ 2 - t53 ^ 2;
t74 = qJD(4) * t53;
t72 = t53 * qJD(3);
t71 = t39 * t82;
t70 = 0.2e1 * qJD(3) * t40;
t69 = t47 * t72;
t68 = t48 * t72;
t67 = t51 * t72;
t33 = t39 * t72;
t65 = t76 * t53;
t64 = t76 * qJD(4);
t63 = 0.2e1 * t67;
t62 = 0.2e1 * t64;
t10 = -t47 * t51 * pkin(7) + t15;
t23 = t48 * t28;
t7 = -pkin(7) * t81 + t23 + (-t39 * t47 - pkin(4)) * t53;
t61 = t52 * t10 + t50 * t7;
t60 = t50 * t10 - t52 * t7;
t17 = t47 * t25;
t13 = -t48 * t66 + t17;
t58 = -t12 * t47 + t13 * t48;
t14 = t23 - t71;
t57 = -t14 * t47 + t15 * t48;
t36 = t77 * t47;
t37 = t77 * t48;
t56 = -t52 * t36 - t50 * t37;
t55 = -t50 * t36 + t52 * t37;
t31 = t52 * t47 + t50 * t48;
t27 = t31 * qJD(5);
t54 = -t53 * t27 - t42 * t85;
t24 = (pkin(4) * t47 + t39) * t51;
t21 = t85 * t51;
t20 = t31 * t51;
t19 = pkin(4) * t69 + t33;
t11 = -t26 * t53 + t31 * t42;
t9 = t51 * t26 + t31 * t72;
t8 = t51 * t27 + t50 * t69 - t52 * t68;
t6 = t17 + (-pkin(7) * t82 - t39 * t81) * qJD(3);
t5 = -t31 * qJD(4) - t55 * qJD(5);
t4 = -qJD(4) * t85 - t56 * qJD(5);
t3 = (pkin(4) * t51 - pkin(7) * t80) * qJD(3) + t12;
t2 = -t61 * qJD(5) + t52 * t3 - t50 * t6;
t1 = t60 * qJD(5) - t50 * t3 - t52 * t6;
t16 = [0, 0, 0, 0, t63, -0.2e1 * t75 * qJD(3), 0, 0, 0, t51 * t70, t53 * t70, -0.2e1 * t12 * t53 + 0.2e1 * (t14 + 0.2e1 * t71) * t42, 0.2e1 * t13 * t53 + 0.2e1 * (-t15 + 0.2e1 * t32) * t42, 0.2e1 * (-t12 * t48 - t13 * t47) * t51 + 0.2e1 * (-t14 * t48 - t15 * t47) * t72, 0.2e1 * t39 ^ 2 * t67 + 0.2e1 * t14 * t12 + 0.2e1 * t15 * t13, -0.2e1 * t21 * t8, 0.2e1 * t8 * t20 - 0.2e1 * t21 * t9, 0.2e1 * t21 * t42 + 0.2e1 * t8 * t53, -0.2e1 * t20 * t42 + 0.2e1 * t53 * t9, -0.2e1 * t67, 0.2e1 * t19 * t20 - 0.2e1 * t2 * t53 + 0.2e1 * t24 * t9 - 0.2e1 * t60 * t42, -0.2e1 * t1 * t53 + 0.2e1 * t19 * t21 - 0.2e1 * t24 * t8 - 0.2e1 * t61 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t51 + (t75 * t39 + t57 * t53) * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t76) * t63, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t72, -t42, 0, -t33, t66, t47 * t74 + (t59 * t47 - t32) * qJD(3), t48 * t74 + (t59 * t48 + t71) * qJD(3), t58, -pkin(3) * t33 + t58 * qJ(4) + t57 * qJD(4), t21 * t26 - t8 * t31, -t26 * t20 - t21 * t27 - t31 * t9 - t8 * t85, t11, -t54, 0, -t19 * t85 + t24 * t27 + t41 * t9 + t56 * t42 - t5 * t53, t19 * t31 + t24 * t26 - t4 * t53 - t41 * t8 - t55 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t72, -t48 * t42, t47 * t42, qJD(3) * t65, t51 * t64 + (qJ(4) * t65 - t83) * qJD(3), 0, 0, 0, 0, 0, t54, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, qJ(4) * t62, 0.2e1 * t31 * t26, 0.2e1 * t26 * t85 - 0.2e1 * t31 * t27, 0, 0, 0, t27 * t84, t26 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t68, 0, t33, 0, 0, 0, 0, 0, t9, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, t42, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t27, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t16;
