% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:25:33
% EndTime: 2019-12-05 18:25:35
% DurationCPUTime: 0.47s
% Computational Cost: add. (513->86), mult. (1263->169), div. (0->0), fcn. (1075->6), ass. (0->74)
t45 = cos(qJ(5));
t40 = t45 ^ 2;
t42 = sin(qJ(5));
t73 = t42 ^ 2 - t40;
t59 = t73 * qJD(5);
t82 = qJD(2) + qJD(4);
t48 = pkin(1) + pkin(2);
t44 = sin(qJ(2));
t74 = pkin(3) + qJ(3);
t29 = t74 * t44;
t47 = cos(qJ(2));
t30 = t74 * t47;
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t16 = -t43 * t29 + t46 * t30;
t60 = qJD(2) * t74;
t50 = t47 * qJD(3) - t44 * t60;
t67 = t44 * qJD(3);
t51 = t47 * t60 + t67;
t7 = t16 * qJD(4) + t43 * t50 + t46 * t51;
t4 = t7 * t42;
t75 = t46 * t29;
t15 = t43 * t30 + t75;
t37 = qJD(5) * t45;
t81 = t15 * t37 + t4;
t26 = t43 * t44 - t46 * t47;
t13 = t82 * t26;
t27 = t43 * t47 + t46 * t44;
t80 = t27 * t13;
t79 = t27 * t45;
t14 = t82 * t27;
t78 = t42 * t14;
t77 = t45 * t13;
t76 = t45 * t14;
t68 = t44 * qJD(2);
t36 = pkin(1) * t68;
t28 = pkin(2) * t68 + t36;
t72 = qJD(4) * t43;
t71 = qJD(4) * t48;
t70 = qJD(5) * t42;
t69 = qJD(5) * t46;
t66 = t47 * qJD(2);
t31 = t48 * t47;
t65 = t43 * t71;
t64 = t46 * t71;
t63 = t42 * t37;
t62 = t44 * t66;
t61 = -0.4e1 * t42 * t79;
t39 = t44 ^ 2;
t41 = t47 ^ 2;
t58 = (t39 + t41) * qJD(3);
t17 = -t27 * pkin(4) - t31;
t57 = t45 * t16 + t42 * t17;
t56 = t42 * t16 - t45 * t17;
t33 = t43 * t48 + pkin(4);
t55 = t27 * t46 * t48 + t26 * t33;
t54 = -t42 * t13 + t27 * t37;
t53 = -t27 * t70 - t77;
t9 = t26 * t37 + t78;
t52 = t26 * t70 - t76;
t49 = -t14 * t33 + (t13 * t46 + (-t26 * t46 + t27 * t43) * qJD(4)) * t48;
t32 = 0.2e1 * t63;
t25 = -0.2e1 * t59;
t24 = t27 ^ 2;
t19 = (-t42 * t69 - t45 * t72) * t48;
t18 = (t42 * t72 - t45 * t69) * t48;
t11 = t15 * t70;
t10 = t13 * pkin(4) + t28;
t6 = qJD(4) * t75 + t30 * t72 + t43 * t51 - t46 * t50;
t5 = -t27 * t59 - t42 * t77;
t3 = qJD(5) * t61 + t73 * t13;
t2 = -t57 * qJD(5) + t45 * t10 + t42 * t6;
t1 = t56 * qJD(5) - t42 * t10 + t45 * t6;
t8 = [0, 0, 0, 0.2e1 * t62, 0.2e1 * (-t39 + t41) * qJD(2), 0, 0, 0, 0, 0, 0.2e1 * t58, -0.2e1 * pkin(1) ^ 2 * t62 + 0.2e1 * qJ(3) * t58, -0.2e1 * t80, 0.2e1 * t13 * t26 - 0.2e1 * t27 * t14, 0, 0, 0, -0.2e1 * t31 * t14 + 0.2e1 * t28 * t26, 0.2e1 * t31 * t13 + 0.2e1 * t28 * t27, -0.2e1 * t24 * t63 - 0.2e1 * t40 * t80, -t13 * t61 + 0.2e1 * t24 * t59, 0.2e1 * t53 * t26 + 0.2e1 * t27 * t76, -0.2e1 * t54 * t26 - 0.2e1 * t27 * t78, 0.2e1 * t26 * t14, -0.2e1 * t56 * t14 + 0.2e1 * t54 * t15 + 0.2e1 * t2 * t26 + 0.2e1 * t27 * t4, 0.2e1 * t1 * t26 - 0.2e1 * t57 * t14 + 0.2e1 * t53 * t15 + 0.2e1 * t7 * t79; 0, 0, 0, 0, 0, t66, -t68, 0, 0, 0, -pkin(1) * t66, (-qJ(3) * t66 - t67) * pkin(1), 0, 0, -t13, -t14, 0, -t7, t6, t5, t3, t9, -t52, 0, t11 + (-t55 * qJD(5) - t7) * t45 + t49 * t42, t49 * t45 + t55 * t70 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t65, -0.2e1 * t64, t32, t25, 0, 0, 0, 0.2e1 * t19, 0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, t14, -t13, 0, 0, 0, 0, 0, -t52, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, -t7, t6, t5, t3, t9, -t52, 0, -t9 * pkin(4) - t7 * t45 + t11, t52 * pkin(4) + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t64, t32, t25, 0, 0, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t54, t14, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t70, 0, -t33 * t37 - t42 * t64, t33 * t70 - t45 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t70, 0, -pkin(4) * t37, pkin(4) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
