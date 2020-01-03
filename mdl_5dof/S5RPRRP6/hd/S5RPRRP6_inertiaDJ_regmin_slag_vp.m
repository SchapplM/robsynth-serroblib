% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP6
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
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:16
% EndTime: 2019-12-31 18:43:18
% DurationCPUTime: 0.58s
% Computational Cost: add. (433->107), mult. (1031->199), div. (0->0), fcn. (749->6), ass. (0->71)
t38 = sin(qJ(3));
t79 = -0.4e1 * t38;
t37 = sin(qJ(4));
t40 = cos(qJ(3));
t62 = t40 * qJD(3);
t39 = cos(qJ(4));
t65 = qJD(4) * t39;
t42 = t37 * t62 + t38 * t65;
t78 = t42 * pkin(4);
t33 = t38 ^ 2;
t48 = (-t40 ^ 2 + t33) * qJD(3);
t32 = t37 ^ 2;
t34 = t39 ^ 2;
t49 = (t32 - t34) * qJD(4);
t29 = -cos(pkin(8)) * pkin(1) - pkin(2);
t46 = -t40 * pkin(3) - t38 * pkin(7);
t20 = t29 + t46;
t45 = pkin(3) * t38 - pkin(7) * t40;
t23 = t45 * qJD(3);
t77 = -t20 * t65 - t37 * t23;
t28 = sin(pkin(8)) * pkin(1) + pkin(6);
t73 = t39 * t40;
t22 = t28 * t73;
t9 = t37 * t20;
t76 = t22 + t9;
t75 = t28 * t37;
t74 = t38 * t39;
t72 = -qJ(5) - pkin(7);
t31 = t38 * qJD(3);
t53 = t28 * t31;
t71 = t39 * t23 + t37 * t53;
t68 = qJ(5) * t38;
t67 = t39 * qJ(5);
t66 = qJD(4) * t37;
t64 = qJD(4) * t40;
t63 = t39 * qJD(5);
t61 = -0.2e1 * pkin(3) * qJD(4);
t60 = 0.2e1 * qJD(3) * t29;
t59 = pkin(4) * t66;
t58 = t37 * t64;
t57 = t39 * t64;
t56 = t32 * t62;
t55 = t37 * t65;
t54 = t38 * t62;
t52 = t39 * t62;
t51 = t28 * t62;
t50 = qJD(4) * t72;
t47 = t37 * t52;
t10 = t39 * t20;
t5 = -t38 * t67 + t10 + (-pkin(4) - t75) * t40;
t6 = -t37 * t68 + t76;
t44 = -t37 * t6 - t39 * t5;
t43 = t37 * t5 - t39 * t6;
t14 = t38 * t66 - t52;
t15 = t39 * t31 + t58;
t11 = t37 * t50 + t63;
t12 = -t37 * qJD(5) + t39 * t50;
t25 = t72 * t37;
t26 = t72 * t39;
t41 = t11 * t39 - t12 * t37 + (-t25 * t39 + t26 * t37) * qJD(4);
t30 = -t39 * pkin(4) - pkin(3);
t27 = t34 * t62;
t24 = t34 * t54;
t17 = t31 * t37 - t57;
t13 = (pkin(4) * t37 + t28) * t38;
t7 = t51 + t78;
t4 = -t76 * qJD(4) + t71;
t3 = t15 * t28 + t77;
t2 = (-qJ(5) * qJD(4) - qJD(3) * t28) * t74 + (-qJD(5) * t38 + (-qJ(5) * qJD(3) - qJD(4) * t28) * t40) * t37 - t77;
t1 = -t38 * t63 + (pkin(4) * t38 - t40 * t67) * qJD(3) + (-t22 + (-t20 + t68) * t37) * qJD(4) + t71;
t8 = [0, 0, 0, 0, 0.2e1 * t54, -0.2e1 * t48, 0, 0, 0, t38 * t60, t40 * t60, -0.2e1 * t33 * t55 + 0.2e1 * t24, 0.2e1 * t33 * t49 + t47 * t79, 0.2e1 * t38 * t58 + 0.2e1 * t39 * t48, -0.2e1 * t37 * t48 + 0.2e1 * t38 * t57, -0.2e1 * t54, 0.2e1 * t10 * t31 - 0.2e1 * t4 * t40 + 0.2e1 * (t33 * t65 + t37 * t54) * t28, -0.2e1 * t33 * t28 * t66 - 0.2e1 * t3 * t40 + 0.2e1 * (t22 - t9) * t31, 0.2e1 * t44 * t62 + 0.2e1 * (qJD(4) * t43 - t1 * t39 - t2 * t37) * t38, 0.2e1 * t5 * t1 + 0.2e1 * t13 * t7 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-qJD(3) * t43 - t7) * t40 + (qJD(3) * t13 + qJD(4) * t44 - t1 * t37 + t2 * t39) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t24 + 0.2e1 * (t32 - 0.1e1) * t54; 0, 0, 0, 0, 0, 0, t62, -t31, 0, -t51, t53, -t38 * t49 + t47, t55 * t79 + t27 - t56, t17, t15, 0, (pkin(7) * t73 + (-pkin(3) * t39 + t75) * t38) * qJD(4) + (t37 * t46 - t22) * qJD(3), (t28 * t74 + t37 * t45) * qJD(4) + (t46 * t39 + t40 * t75) * qJD(3), (-t25 * t62 - t12 * t38 + t2 + (t26 * t38 - t5) * qJD(4)) * t39 + (t26 * t62 - t11 * t38 - t1 + (t25 * t38 - t6) * qJD(4)) * t37, t1 * t25 + t6 * t11 + t5 * t12 + t13 * t59 - t2 * t26 + t7 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t62, 0, 0, 0, 0, 0, -t15, t17, t27 + t56, (-t59 + (-t25 * t37 - t26 * t39) * qJD(3)) * t40 + (qJD(3) * t30 + t41) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t55, -0.2e1 * t49, 0, 0, 0, t37 * t61, t39 * t61, 0.2e1 * t41, -0.2e1 * t26 * t11 + 0.2e1 * t25 * t12 + 0.2e1 * t30 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t42, t31, t4, t3, t14 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t14, 0, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t66, 0, -pkin(7) * t65, pkin(7) * t66, -pkin(4) * t65, t12 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
