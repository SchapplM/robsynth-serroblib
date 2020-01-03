% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP4
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
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:24
% EndTime: 2020-01-03 11:50:26
% DurationCPUTime: 0.49s
% Computational Cost: add. (683->93), mult. (1686->177), div. (0->0), fcn. (1486->6), ass. (0->66)
t45 = sin(qJ(4));
t46 = sin(qJ(3));
t47 = cos(qJ(4));
t48 = cos(qJ(3));
t28 = t45 * t48 + t47 * t46;
t43 = sin(pkin(8));
t22 = t28 * t43;
t64 = qJD(3) + qJD(4);
t44 = cos(pkin(8));
t29 = -t44 * pkin(2) - t43 * pkin(6) - pkin(1);
t79 = pkin(7) * t43;
t53 = -t29 + t79;
t58 = t48 * t44 * qJ(2);
t82 = t53 * t46 - t58;
t81 = 0.2e1 * t44;
t80 = 2 * qJD(3);
t78 = t47 * pkin(3);
t77 = t45 * t82;
t76 = t45 * t46;
t75 = t47 * t82;
t74 = t47 * t48;
t68 = qJD(3) * t48;
t70 = qJD(2) * t48;
t73 = t29 * t68 + t44 * t70;
t56 = t43 * t68;
t25 = pkin(3) * t56 + t43 * qJD(2);
t26 = (pkin(3) * t46 + qJ(2)) * t43;
t72 = qJ(2) * t46;
t71 = qJD(2) * t46;
t69 = qJD(3) * t46;
t67 = qJD(4) * t45;
t66 = qJD(4) * t47;
t65 = qJ(2) * qJD(3);
t13 = (-t44 * t72 - t48 * t79) * qJD(3) + t73;
t55 = t44 * t71;
t14 = qJD(3) * t82 - t55;
t16 = -t53 * t48 + (-pkin(3) - t72) * t44;
t63 = -t47 * t13 - t45 * t14 - t16 * t66;
t62 = t43 * t76;
t61 = t43 * t74;
t60 = pkin(3) * t67;
t59 = pkin(3) * t66;
t57 = t43 * t69;
t54 = t46 * t65;
t52 = -t45 * t13 + t47 * t14;
t51 = t43 * t44 * t80;
t41 = t43 ^ 2;
t50 = 0.2e1 * (t44 ^ 2 + t41) * qJD(2);
t49 = -t45 * t16 + t75;
t3 = -t67 * t82 + t63;
t4 = t49 * qJD(4) + t52;
t38 = pkin(4) + t78;
t27 = t74 - t76;
t23 = t61 - t62;
t20 = t64 * t28;
t19 = -t47 * t68 - t48 * t66 + t64 * t76;
t18 = -t55 + (-t46 * t29 - t58) * qJD(3);
t17 = t44 * t54 - t73;
t12 = -qJD(4) * t62 - t45 * t57 + t64 * t61;
t11 = t64 * t22;
t7 = t12 * pkin(4) + t25;
t6 = -t22 * qJ(5) - t49;
t5 = -t44 * pkin(4) - t23 * qJ(5) + t47 * t16 + t77;
t2 = t11 * qJ(5) - t23 * qJD(5) + t4;
t1 = -t12 * qJ(5) - t22 * qJD(5) - t3;
t8 = [0, 0, 0, 0, 0, t50, qJ(2) * t50, -0.2e1 * t41 * t46 * t68, (t46 ^ 2 - t48 ^ 2) * t41 * t80, t46 * t51, t48 * t51, 0, -0.2e1 * t18 * t44 + 0.2e1 * (t48 * t65 + t71) * t41, -0.2e1 * t17 * t44 + 0.2e1 * (-t54 + t70) * t41, -0.2e1 * t23 * t11, 0.2e1 * t11 * t22 - 0.2e1 * t23 * t12, t11 * t81, t12 * t81, 0, 0.2e1 * t26 * t12 + 0.2e1 * t25 * t22 - 0.2e1 * t4 * t44, -0.2e1 * t26 * t11 + 0.2e1 * t25 * t23 - 0.2e1 * t3 * t44, -0.2e1 * t1 * t22 + 0.2e1 * t5 * t11 - 0.2e1 * t6 * t12 - 0.2e1 * t2 * t23, 0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * (t22 * pkin(4) + t26) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t69, t44 * t68, 0, 0, 0, 0, 0, t20 * t44, -t19 * t44, t27 * t11 - t28 * t12 + t19 * t22 + t20 * t23, t1 * t28 - t6 * t19 + t2 * t27 - t5 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t28 * t19 - 0.2e1 * t27 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t56, 0, t18, t17, 0, 0, -t11, -t12, 0, (t75 + (t44 * pkin(3) - t16) * t45) * qJD(4) + t52, (t44 * t78 - t77) * qJD(4) + t63, t38 * t11 + (-t12 * t45 + (-t22 * t47 + t23 * t45) * qJD(4)) * pkin(3), t2 * t38 + (t1 * t45 + (-t45 * t5 + t47 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t68, 0, 0, 0, 0, 0, -t20, t19, 0, -t20 * t38 + (-t19 * t45 + (-t27 * t45 + t28 * t47) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t60, -0.2e1 * t59, 0, 0.2e1 * (-t38 + t78) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, t4, t3, pkin(4) * t11, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t19, 0, -t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t59, 0, -pkin(4) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
