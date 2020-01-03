% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:36
% EndTime: 2019-12-31 18:45:39
% DurationCPUTime: 0.67s
% Computational Cost: add. (555->123), mult. (1277->209), div. (0->0), fcn. (925->6), ass. (0->78)
t37 = sin(qJ(3));
t38 = cos(qJ(4));
t36 = sin(qJ(4));
t46 = pkin(4) * t38 + qJ(5) * t36;
t88 = qJD(4) * t46 - t38 * qJD(5);
t91 = t88 * t37;
t90 = -0.4e1 * t37;
t29 = -cos(pkin(8)) * pkin(1) - pkin(2);
t39 = cos(qJ(3));
t85 = t37 * pkin(7);
t51 = -pkin(3) * t39 - t85;
t18 = t29 + t51;
t10 = t36 * t18;
t28 = sin(pkin(8)) * pkin(1) + pkin(6);
t80 = t28 * t39;
t20 = t38 * t80;
t43 = t20 + t10;
t31 = t36 ^ 2;
t33 = t38 ^ 2;
t54 = (t31 - t33) * qJD(4);
t45 = pkin(4) * t36 - qJ(5) * t38;
t11 = qJD(4) * t45 - t36 * qJD(5);
t22 = -pkin(3) - t46;
t42 = t28 + t45;
t8 = t42 * t37;
t89 = (-t22 * t39 + t85) * qJD(3) - qJD(4) * t8 - t11 * t37;
t86 = pkin(7) * t39;
t50 = pkin(3) * t37 - t86;
t21 = t50 * qJD(3);
t71 = t37 * qJD(3);
t57 = t28 * t71;
t4 = -qJD(4) * t43 + t38 * t21 + t36 * t57;
t87 = 0.2e1 * qJD(5);
t30 = qJD(4) * t38;
t84 = t18 * t30 + t36 * t21;
t83 = t22 * t37;
t82 = t28 * t36;
t81 = t28 * t38;
t32 = t37 ^ 2;
t77 = -t39 ^ 2 + t32;
t76 = qJD(3) * t8;
t75 = qJD(3) * t38;
t74 = qJD(4) * t36;
t73 = qJD(4) * t37;
t72 = qJD(4) * t39;
t69 = t39 * qJD(3);
t68 = -0.2e1 * pkin(3) * qJD(4);
t67 = 0.2e1 * qJD(3) * t29;
t66 = pkin(4) * t71;
t65 = pkin(7) * t74;
t64 = pkin(7) * t30;
t63 = t28 * t74;
t62 = t36 * t72;
t61 = t38 * t72;
t60 = t31 * t69;
t59 = t36 * t30;
t58 = t37 * t69;
t56 = t38 * t71;
t55 = t38 * t69;
t53 = t77 * qJD(3);
t52 = t36 * t55;
t6 = -qJ(5) * t39 + t43;
t7 = -t38 * t18 + (pkin(4) + t82) * t39;
t49 = t36 * t7 + t38 * t6;
t48 = -t36 * t6 + t38 * t7;
t13 = t56 + t62;
t1 = (-qJD(5) - t63) * t39 + (qJ(5) - t81) * t71 + t84;
t2 = -t4 - t66;
t40 = qJD(4) * t48 + t1 * t38 + t2 * t36;
t26 = t33 * t69;
t24 = pkin(7) * t61;
t23 = t33 * t58;
t15 = -t36 * t71 + t61;
t14 = t30 * t37 + t36 * t69;
t12 = t36 * t73 - t55;
t5 = t42 * t69 + t91;
t3 = t13 * t28 - t84;
t9 = [0, 0, 0, 0, 0.2e1 * t58, -0.2e1 * t53, 0, 0, 0, t37 * t67, t39 * t67, -0.2e1 * t32 * t59 + 0.2e1 * t23, 0.2e1 * t32 * t54 + t52 * t90, 0.2e1 * t37 * t62 + 0.2e1 * t75 * t77, -0.2e1 * t36 * t53 + 0.2e1 * t37 * t61, -0.2e1 * t58, 0.2e1 * t18 * t56 - 0.2e1 * t4 * t39 + 0.2e1 * (t30 * t32 + t36 * t58) * t28, -0.2e1 * t32 * t63 - 0.2e1 * t3 * t39 + 0.2e1 * (-t10 + t20) * t71, 0.2e1 * (t36 * t76 + t2) * t39 + 0.2e1 * (-qJD(3) * t7 + t30 * t8 + t5 * t36) * t37, 0.2e1 * t48 * t69 + 0.2e1 * (-qJD(4) * t49 - t1 * t36 + t2 * t38) * t37, 0.2e1 * (-t75 * t8 - t1) * t39 + 0.2e1 * (qJD(3) * t6 - t5 * t38 + t74 * t8) * t37, 0.2e1 * t1 * t6 + 0.2e1 * t2 * t7 + 0.2e1 * t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (qJD(3) * t49 - t5) * t39 + (t40 + t76) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t23 + 0.2e1 * (t31 - 0.1e1) * t58; 0, 0, 0, 0, 0, 0, t69, -t71, 0, -t28 * t69, t57, -t37 * t54 + t52, t59 * t90 + t26 - t60, -t15, t13, 0, t24 + (-pkin(3) * t38 + t82) * t73 + (t36 * t51 - t20) * qJD(3), (t36 * t50 + t37 * t81) * qJD(4) + (t36 * t80 + t38 * t51) * qJD(3), t24 + (t22 * t73 - t5) * t38 - t89 * t36, t40, (-t5 + (t83 + t86) * qJD(4)) * t36 + t89 * t38, pkin(7) * t40 + t8 * t11 + t5 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t69, 0, 0, 0, 0, 0, -t13, -t15, -t13, t26 + t60, t15, -t39 * t11 + (t83 + (t31 + t33) * t86) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t59, -0.2e1 * t54, 0, 0, 0, t36 * t68, t38 * t68, -0.2e1 * t11 * t38 + 0.2e1 * t22 * t74, 0, -0.2e1 * t11 * t36 - 0.2e1 * t22 * t30, 0.2e1 * t22 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t14, t71, t4, t3, t4 + 0.2e1 * t66, (-pkin(4) * t69 - qJ(5) * t73) * t38 + (-qJ(5) * t69 + (pkin(4) * qJD(4) - qJD(5)) * t37) * t36, (-0.2e1 * qJD(5) - t63) * t39 + (0.2e1 * qJ(5) - t81) * t71 + t84, -pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t12, -t14, 0, -t12, -t45 * t69 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t74, 0, -t64, t65, -t64, -t88, -t65, -t88 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, qJ(5) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t12, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
