% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPPRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:42
% EndTime: 2019-03-09 01:37:43
% DurationCPUTime: 0.52s
% Computational Cost: add. (276->83), mult. (603->172), div. (0->0), fcn. (479->6), ass. (0->68)
t31 = sin(qJ(5));
t23 = t31 ^ 2;
t33 = cos(qJ(5));
t49 = (-t33 ^ 2 + t23) * qJD(5);
t32 = cos(qJ(6));
t24 = t32 ^ 2;
t30 = sin(qJ(6));
t72 = t30 ^ 2 - t24;
t50 = t72 * qJD(6);
t26 = sin(pkin(9));
t27 = cos(pkin(9));
t16 = -qJD(2) * t26 + qJD(3) * t27;
t76 = t23 * t16;
t75 = t30 * t33;
t74 = t32 * t33;
t28 = pkin(3) + qJ(2);
t29 = -pkin(1) - qJ(3);
t73 = t26 * t28 + t27 * t29;
t70 = qJD(6) * t30;
t69 = qJD(6) * t31;
t68 = qJD(6) * t32;
t67 = qJD(6) * t33;
t66 = t31 * qJD(5);
t65 = t33 * qJD(5);
t64 = qJ(2) * qJD(2);
t63 = -0.2e1 * pkin(5) * qJD(6);
t62 = pkin(8) * t68;
t61 = t30 * t67;
t60 = t32 * t67;
t59 = t30 * t66;
t58 = t30 * t68;
t57 = t31 * t65;
t17 = t26 * t66;
t56 = t27 * t66;
t55 = t32 * t66;
t54 = t32 * t65;
t53 = t26 * t65;
t52 = t27 * t65;
t51 = -t26 * t29 + t27 * t28;
t48 = t30 * t17;
t47 = t30 * t56;
t46 = t26 * t55;
t45 = t27 * t55;
t44 = t31 * t54;
t13 = -pkin(4) - t51;
t14 = pkin(7) + t73;
t15 = t27 * qJD(2) + t26 * qJD(3);
t41 = pkin(5) * t31 - pkin(8) * t33;
t43 = -qJD(5) * t41 + t14 * t67 + t15;
t42 = -pkin(5) * t33 - pkin(8) * t31;
t40 = -t14 * t65 + t31 * t16;
t39 = t14 * t66 + t16 * t33;
t9 = t30 * t69 - t54;
t38 = t30 * t65 + t31 * t68;
t37 = -t23 * t70 + t44;
t36 = t23 * t68 + t30 * t57;
t7 = t13 + t42;
t35 = -qJD(6) * t7 + t39;
t34 = 0.2e1 * qJD(2);
t12 = t59 - t60;
t10 = t55 + t61;
t6 = t47 + (-t26 * t30 - t27 * t74) * qJD(6);
t5 = t48 + (-t26 * t74 + t27 * t30) * qJD(6);
t4 = t45 + (-t26 * t32 + t27 * t75) * qJD(6);
t3 = t46 + (t26 * t75 + t27 * t32) * qJD(6);
t2 = t30 * t35 - t32 * t43;
t1 = t30 * t43 + t32 * t35;
t8 = [0, 0, 0, 0, t34, 0.2e1 * t64, t34, 0.2e1 * qJD(3), -0.2e1 * qJD(3) * t29 + 0.2e1 * t64, 0.2e1 * t15, 0.2e1 * t16, 0.2e1 * t15 * t51 - 0.2e1 * t16 * t73, 0.2e1 * t57, -0.2e1 * t49, 0, 0, 0, 0.2e1 * t13 * t66 + 0.2e1 * t15 * t33, 0.2e1 * t13 * t65 - 0.2e1 * t15 * t31, -0.2e1 * t23 * t58 + 0.2e1 * t24 * t57, 0.2e1 * t23 * t50 - 0.4e1 * t30 * t44, 0.2e1 * t31 * t61 + 0.2e1 * t32 * t49, -0.2e1 * t30 * t49 + 0.2e1 * t31 * t60, -0.2e1 * t57, 0.2e1 * t14 * t36 - 0.2e1 * t2 * t33 - 0.2e1 * t30 * t76 + 0.2e1 * t55 * t7, -0.2e1 * t1 * t33 + 0.2e1 * t14 * t37 - 0.2e1 * t32 * t76 - 0.2e1 * t59 * t7; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, -t15 * t26 - t16 * t27, 0, 0, 0, 0, 0, t17, t53, 0, 0, 0, 0, 0, t27 * t36 - t6 * t33 + t46, t27 * t37 - t4 * t33 - t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, t15 * t27 - t16 * t26, 0, 0, 0, 0, 0, -t56, -t52, 0, 0, 0, 0, 0, t26 * t36 - t5 * t33 - t45, t26 * t37 - t3 * t33 + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t66, 0, t40, t39, t30 * t54 - t31 * t50, -0.4e1 * t31 * t58 - t65 * t72, t12, t10, 0 (t62 + (-pkin(5) * t30 - t14 * t32) * qJD(5)) * t33 + (-pkin(8) * qJD(5) * t30 + t16 * t32 + (-pkin(5) * t32 + t14 * t30) * qJD(6)) * t31 (qJD(5) * t42 + t14 * t69) * t32 + (qJD(6) * t41 - t40) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t56, 0, 0, 0, 0, 0, t9 * t27, t38 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t17, 0, 0, 0, 0, 0, t9 * t26, t38 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, 0, 0, 0, 0, -t10, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t58, -0.2e1 * t50, 0, 0, 0, t30 * t63, t32 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t38, t66, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t70, 0, -t62, pkin(8) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
