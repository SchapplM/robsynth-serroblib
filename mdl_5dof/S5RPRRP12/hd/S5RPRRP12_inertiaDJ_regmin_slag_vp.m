% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:20
% EndTime: 2019-12-31 18:57:22
% DurationCPUTime: 0.60s
% Computational Cost: add. (421->109), mult. (954->201), div. (0->0), fcn. (663->4), ass. (0->74)
t32 = sin(qJ(3));
t34 = cos(qJ(3));
t77 = t34 * pkin(7);
t41 = t32 * pkin(3) - t77;
t21 = qJ(2) + t41;
t31 = sin(qJ(4));
t17 = t31 * t21;
t33 = cos(qJ(4));
t35 = -pkin(1) - pkin(6);
t25 = t33 * t32 * t35;
t79 = -t17 - t25;
t28 = t32 ^ 2;
t30 = t34 ^ 2;
t46 = (t28 - t30) * qJD(3);
t27 = t31 ^ 2;
t29 = t33 ^ 2;
t73 = t27 - t29;
t47 = t73 * qJD(4);
t78 = 2 * qJD(2);
t76 = t31 * t35;
t75 = t34 * t35;
t74 = -qJ(5) - pkin(7);
t72 = t27 + t29;
t70 = t28 + t30;
t69 = qJ(5) * t34;
t68 = qJD(4) * t31;
t67 = qJD(4) * t33;
t66 = qJD(4) * t34;
t65 = qJD(4) * t35;
t64 = t32 * qJD(3);
t63 = t33 * qJD(5);
t62 = t34 * qJD(3);
t61 = qJ(2) * qJD(3);
t60 = -0.2e1 * pkin(3) * qJD(4);
t42 = pkin(3) * t34 + pkin(7) * t32;
t19 = t42 * qJD(3) + qJD(2);
t50 = t35 * t62;
t59 = t31 * t19 + t21 * t67 + t33 * t50;
t58 = pkin(4) * t68;
t57 = t31 * t66;
t56 = t31 * t65;
t55 = t33 * t66;
t54 = t31 * t67;
t53 = t35 * t64;
t52 = t33 * t64;
t51 = t32 * t62;
t49 = pkin(4) - t76;
t48 = qJD(4) * t74;
t45 = 0.2e1 * t51;
t44 = t31 * t50;
t43 = t31 * t52;
t18 = t33 * t21;
t5 = t49 * t32 - t33 * t69 + t18;
t6 = -t31 * t69 - t79;
t40 = t31 * t6 + t33 * t5;
t39 = t31 * t5 - t33 * t6;
t38 = t52 + t57;
t14 = t31 * t62 + t32 * t67;
t15 = t31 * t64 - t55;
t37 = t79 * qJD(4) + t33 * t19;
t23 = t74 * t31;
t24 = t74 * t33;
t8 = t31 * t48 + t63;
t9 = -t31 * qJD(5) + t33 * t48;
t36 = -t9 * t31 + t8 * t33 + (-t23 * t33 + t24 * t31) * qJD(4);
t26 = -t33 * pkin(4) - pkin(3);
t20 = (pkin(4) * t31 - t35) * t34;
t12 = t32 * t68 - t33 * t62;
t7 = -t15 * pkin(4) + t53;
t4 = t37 - t44;
t3 = t32 * t56 - t59;
t2 = -qJ(5) * t55 + (-qJD(5) * t34 + (qJ(5) * qJD(3) - t65) * t32) * t31 + t59;
t1 = qJ(5) * t52 + (qJ(5) * t68 + t49 * qJD(3) - t63) * t34 + t37;
t10 = [0, 0, 0, 0, t78, qJ(2) * t78, -0.2e1 * t51, 0.2e1 * t46, 0, 0, 0, 0.2e1 * qJD(2) * t32 + 0.2e1 * t34 * t61, 0.2e1 * qJD(2) * t34 - 0.2e1 * t32 * t61, -0.2e1 * t29 * t51 - 0.2e1 * t30 * t54, 0.2e1 * t30 * t47 + 0.4e1 * t34 * t43, -0.2e1 * t32 * t57 - 0.2e1 * t33 * t46, 0.2e1 * t31 * t46 - 0.2e1 * t32 * t55, t45, -0.2e1 * t30 * t33 * t65 + 0.2e1 * t18 * t62 + 0.2e1 * (t4 + t44) * t32, 0.2e1 * t30 * t56 + 0.2e1 * t3 * t32 + 0.2e1 * (-t17 + t25) * t62, 0.2e1 * t40 * t64 + 0.2e1 * (t39 * qJD(4) - t1 * t33 - t2 * t31) * t34, 0.2e1 * t5 * t1 + 0.2e1 * t6 * t2 + 0.2e1 * t20 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70 * t67, t70 * t68, 0, (-t39 * qJD(3) - t7) * t34 + (qJD(3) * t20 - t40 * qJD(4) - t1 * t31 + t2 * t33) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t72) * t45; 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t62, 0, -t53, -t50, -t34 * t47 - t43, -0.4e1 * t34 * t54 + t73 * t64, t14, -t12, 0, (-t31 * t75 - t42 * t33) * qJD(4) + (t41 * t31 - t25) * qJD(3), (t42 * t31 - t33 * t75) * qJD(4) + (-t33 * t77 + (pkin(3) * t33 + t76) * t32) * qJD(3), (t23 * t64 - t34 * t9 + t2 + (t24 * t34 - t5) * qJD(4)) * t33 + (-t24 * t64 - t34 * t8 - t1 + (t23 * t34 - t6) * qJD(4)) * t31, t1 * t23 - t2 * t24 + t20 * t58 + t7 * t26 + t5 * t9 + t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t62, 0, 0, 0, 0, 0, -t38, t15, t72 * t62, (-t58 + (-t23 * t31 - t24 * t33) * qJD(3)) * t34 + (qJD(3) * t26 + t36) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t54, -0.2e1 * t47, 0, 0, 0, t31 * t60, t33 * t60, 0.2e1 * t36, 0.2e1 * t23 * t9 - 0.2e1 * t24 * t8 + 0.2e1 * t26 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t15, t62, t4, t3, t38 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t12, 0, -t14 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t68, 0, -pkin(7) * t67, pkin(7) * t68, -pkin(4) * t67, t9 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
