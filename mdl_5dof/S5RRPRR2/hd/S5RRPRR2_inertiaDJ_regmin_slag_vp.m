% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:35
% EndTime: 2021-01-15 21:23:39
% DurationCPUTime: 0.52s
% Computational Cost: add. (1365->89), mult. (3068->175), div. (0->0), fcn. (3027->8), ass. (0->70)
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t59 = cos(pkin(9));
t73 = t59 * pkin(2) + pkin(3);
t58 = sin(pkin(9));
t84 = pkin(2) * t58;
t85 = t61 * t84 - t64 * t73;
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t81 = -qJ(3) - pkin(6);
t70 = qJD(2) * t81;
t41 = t65 * qJD(3) + t62 * t70;
t42 = -t62 * qJD(3) + t65 * t70;
t23 = -t58 * t41 + t59 * t42;
t78 = t65 * qJD(2);
t79 = t62 * qJD(2);
t45 = -t58 * t79 + t59 * t78;
t20 = -t45 * pkin(7) + t23;
t24 = t59 * t41 + t58 * t42;
t47 = t58 * t65 + t59 * t62;
t44 = t47 * qJD(2);
t21 = -t44 * pkin(7) + t24;
t52 = t81 * t62;
t53 = t81 * t65;
t29 = t59 * t52 + t58 * t53;
t25 = -t47 * pkin(7) + t29;
t30 = t58 * t52 - t59 * t53;
t46 = t58 * t62 - t59 * t65;
t26 = -t46 * pkin(7) + t30;
t67 = t64 * t25 - t61 * t26;
t7 = -t67 * qJD(4) - t61 * t20 - t64 * t21;
t39 = pkin(4) - t85;
t83 = -pkin(4) - t39;
t40 = t61 * t73 + t64 * t84;
t63 = cos(qJ(5));
t82 = t40 * t63;
t60 = sin(qJ(5));
t80 = qJD(5) * t60;
t76 = -0.2e1 * pkin(1) * qJD(2);
t57 = pkin(2) * t79;
t75 = pkin(4) * t80;
t74 = qJD(5) * t63 * pkin(4);
t56 = -t65 * pkin(2) - pkin(1);
t31 = t44 * pkin(3) + t57;
t36 = t85 * qJD(4);
t37 = t40 * qJD(4);
t72 = t60 * t36 - t63 * t37;
t71 = t60 * t37 + t40 * t80;
t66 = -t61 * t25 - t64 * t26;
t27 = t64 * t46 + t61 * t47;
t28 = -t61 * t46 + t64 * t47;
t17 = t63 * t27 + t60 * t28;
t18 = -t60 * t27 + t63 * t28;
t34 = t46 * pkin(3) + t56;
t8 = t66 * qJD(4) + t64 * t20 - t61 * t21;
t22 = t27 * pkin(4) + t34;
t16 = t28 * qJD(4) + t64 * t44 + t61 * t45;
t15 = -t27 * qJD(4) - t61 * t44 + t64 * t45;
t13 = (-t39 * t60 - t82) * qJD(5) + t72;
t12 = (-qJD(5) * t39 + t36) * t63 + t71;
t11 = t16 * pkin(4) + t31;
t10 = -t27 * pkin(8) - t66;
t9 = -t28 * pkin(8) + t67;
t6 = t18 * qJD(5) + t60 * t15 + t63 * t16;
t5 = -t17 * qJD(5) + t63 * t15 - t60 * t16;
t4 = -t15 * pkin(8) + t8;
t3 = -t16 * pkin(8) - t7;
t2 = -t60 * t3 + t63 * t4 + (-t10 * t63 - t60 * t9) * qJD(5);
t1 = -t63 * t3 - t60 * t4 + (t10 * t60 - t63 * t9) * qJD(5);
t14 = [0, 0, 0, 0.2e1 * t62 * t78, 0.2e1 * (-t62 ^ 2 + t65 ^ 2) * qJD(2), 0, 0, 0, t62 * t76, t65 * t76, 0.2e1 * t56 * t44 + 0.2e1 * t46 * t57, 0.2e1 * t56 * t45 + 0.2e1 * t47 * t57, -0.2e1 * t23 * t47 - 0.2e1 * t24 * t46 - 0.2e1 * t29 * t45 - 0.2e1 * t30 * t44, 0.2e1 * t29 * t23 + 0.2e1 * t30 * t24 + 0.2e1 * t56 * t57, 0.2e1 * t28 * t15, -0.2e1 * t15 * t27 - 0.2e1 * t28 * t16, 0, 0, 0, 0.2e1 * t34 * t16 + 0.2e1 * t31 * t27, 0.2e1 * t34 * t15 + 0.2e1 * t31 * t28, 0.2e1 * t18 * t5, -0.2e1 * t5 * t17 - 0.2e1 * t18 * t6, 0, 0, 0, 0.2e1 * t11 * t17 + 0.2e1 * t22 * t6, 0.2e1 * t11 * t18 + 0.2e1 * t22 * t5; 0, 0, 0, 0, 0, t78, -t79, 0, -pkin(6) * t78, pkin(6) * t79, t23, -t24, (-t44 * t58 - t45 * t59) * pkin(2), (t23 * t59 + t24 * t58) * pkin(2), 0, 0, t15, -t16, 0, t8, t7, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t37, 0.2e1 * t36, 0, 0, 0, 0, 0, 0.2e1 * t13, 0.2e1 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, 0, t57, 0, 0, 0, 0, 0, t16, t15, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, t8, t7, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t36, 0, 0, 0, 0, 0, (t83 * t60 - t82) * qJD(5) + t72, (t83 * qJD(5) + t36) * t63 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t75, -0.2e1 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
