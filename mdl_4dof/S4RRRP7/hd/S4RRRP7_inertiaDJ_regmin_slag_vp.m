% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:04
% EndTime: 2019-12-31 17:21:06
% DurationCPUTime: 0.52s
% Computational Cost: add. (356->93), mult. (960->197), div. (0->0), fcn. (648->4), ass. (0->67)
t29 = cos(qJ(2));
t27 = sin(qJ(2));
t75 = t27 * pkin(6);
t42 = -t29 * pkin(2) - t75;
t17 = -pkin(1) + t42;
t28 = cos(qJ(3));
t77 = pkin(5) * t29;
t20 = t28 * t77;
t26 = sin(qJ(3));
t82 = t26 * t17 + t20;
t24 = t28 ^ 2;
t71 = t26 ^ 2 - t24;
t46 = t71 * qJD(3);
t76 = pkin(6) * t29;
t41 = pkin(2) * t27 - t76;
t15 = t41 * qJD(2);
t21 = qJD(3) * t28;
t65 = t27 * qJD(2);
t49 = t28 * t65;
t66 = qJD(3) * t29;
t55 = t26 * t66;
t33 = t49 + t55;
t4 = pkin(5) * t33 - t26 * t15 - t17 * t21;
t38 = t28 * pkin(3) + t26 * qJ(4);
t16 = -pkin(2) - t38;
t37 = pkin(3) * t26 - qJ(4) * t28;
t34 = pkin(5) + t37;
t8 = t34 * t27;
t9 = qJD(3) * t37 - t26 * qJD(4);
t81 = qJD(2) * (-t16 * t29 + t75) - qJD(3) * t8 - t27 * t9;
t53 = t26 * t65;
t5 = pkin(5) * t53 - qJD(3) * t82 + t28 * t15;
t80 = qJD(3) * t38 - t28 * qJD(4);
t79 = 0.2e1 * qJD(4);
t78 = pkin(5) * t26;
t23 = t27 ^ 2;
t70 = -t29 ^ 2 + t23;
t69 = qJD(2) * t28;
t68 = qJD(3) * t26;
t67 = qJD(3) * t27;
t63 = t29 * qJD(2);
t62 = t29 * qJD(4);
t61 = qJ(4) * qJD(2);
t60 = -0.2e1 * pkin(1) * qJD(2);
t59 = -0.2e1 * pkin(2) * qJD(3);
t58 = pkin(3) * t65;
t57 = pkin(6) * t68;
t56 = pkin(6) * t21;
t54 = t28 * t66;
t52 = t26 * t21;
t51 = t27 * t63;
t50 = t27 * t21;
t48 = t28 * t63;
t47 = t27 * t61;
t45 = t70 * qJD(2);
t44 = 0.2e1 * t51;
t43 = t26 * t48;
t6 = -t29 * qJ(4) + t82;
t7 = -t28 * t17 + (pkin(3) + t78) * t29;
t40 = -t26 * t6 + t28 * t7;
t1 = -t4 + t47 - t62;
t2 = -t5 - t58;
t30 = qJD(3) * t40 + t1 * t28 + t2 * t26;
t19 = pkin(6) * t54;
t10 = -t26 * t67 + t48;
t3 = t80 * t27 + t34 * t63;
t11 = [0, 0, 0, t44, -0.2e1 * t45, 0, 0, 0, t27 * t60, t29 * t60, -0.2e1 * t23 * t52 + 0.2e1 * t24 * t51, 0.2e1 * t23 * t46 - 0.4e1 * t27 * t43, 0.2e1 * t27 * t55 + 0.2e1 * t69 * t70, -0.2e1 * t26 * t45 + 0.2e1 * t29 * t50, -0.2e1 * t51, 0.2e1 * t17 * t49 - 0.2e1 * t5 * t29 + 0.2e1 * (t21 * t23 + t26 * t51) * pkin(5), -0.2e1 * t4 * t29 - 0.2e1 * t82 * t65 + 0.2e1 * (-t23 * t68 + t28 * t44) * pkin(5), 0.2e1 * (qJD(2) * t26 * t8 + t2) * t29 + 0.2e1 * (-qJD(2) * t7 + t21 * t8 + t3 * t26) * t27, 0.2e1 * t40 * t63 + 0.2e1 * (-t1 * t26 + t2 * t28 + (-t26 * t7 - t28 * t6) * qJD(3)) * t27, 0.2e1 * (-t69 * t8 - t1) * t29 + 0.2e1 * (qJD(2) * t6 - t3 * t28 + t68 * t8) * t27, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, t63, -t65, 0, -pkin(5) * t63, pkin(5) * t65, -t27 * t46 + t43, -0.4e1 * t26 * t50 - t63 * t71, t53 - t54, t33, 0, t19 + (-pkin(2) * t28 + t78) * t67 + (t26 * t42 - t20) * qJD(2), (pkin(5) * t27 * t28 + t26 * t41) * qJD(3) + (t26 * t77 + t42 * t28) * qJD(2), t19 + (t16 * t67 - t3) * t28 - t81 * t26, t30, (-t3 + (t16 * t27 + t76) * qJD(3)) * t26 + t81 * t28, pkin(6) * t30 + t3 * t16 + t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t52, -0.2e1 * t46, 0, 0, 0, t26 * t59, t28 * t59, 0.2e1 * t16 * t68 - 0.2e1 * t9 * t28, 0, -0.2e1 * t16 * t21 - 0.2e1 * t9 * t26, 0.2e1 * t16 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t26 * t63 - t50, t65, t5, t4, t5 + 0.2e1 * t58, (-pkin(3) * t63 - qJ(4) * t67) * t28 + (-t29 * t61 + (pkin(3) * qJD(3) - qJD(4)) * t27) * t26, -t4 + 0.2e1 * t47 - 0.2e1 * t62, -t2 * pkin(3) + t1 * qJ(4) + t6 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t68, 0, -t56, t57, -t56, -t80, -t57, -t80 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, qJ(4) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t10, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
