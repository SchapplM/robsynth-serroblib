% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:34
% EndTime: 2019-12-05 16:37:38
% DurationCPUTime: 0.88s
% Computational Cost: add. (576->152), mult. (1751->314), div. (0->0), fcn. (1585->10), ass. (0->91)
t49 = sin(qJ(5));
t52 = cos(qJ(5));
t79 = qJ(4) * qJD(5);
t102 = -qJD(4) * t52 + t49 * t79;
t45 = sin(pkin(10));
t101 = 0.2e1 * t45;
t100 = 2 * qJD(5);
t48 = cos(pkin(5));
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t46 = sin(pkin(5));
t51 = sin(qJ(2));
t97 = t46 * t51;
t29 = t48 * t50 + t53 * t97;
t54 = cos(qJ(2));
t96 = t46 * t54;
t74 = qJD(2) * t96;
t18 = t29 * qJD(3) + t50 * t74;
t99 = t18 * t50;
t98 = t45 * t53;
t27 = -t50 * qJD(4) + (pkin(3) * t50 - qJ(4) * t53) * qJD(3);
t47 = cos(pkin(10));
t95 = t47 * t27;
t63 = -t53 * pkin(3) - t50 * qJ(4);
t36 = -pkin(2) + t63;
t94 = t47 * t36;
t93 = t47 * t50;
t92 = t47 * t53;
t91 = t49 * t50;
t90 = t52 * t53;
t41 = pkin(7) * t92;
t24 = t45 * t36 + t41;
t89 = qJ(4) * t47;
t88 = qJD(2) * t51;
t87 = qJD(4) * t49;
t85 = qJD(4) * t53;
t84 = qJD(5) * t49;
t83 = qJD(5) * t52;
t71 = pkin(7) * t45 + pkin(4);
t20 = t71 * t53 - t94;
t82 = t20 * qJD(5);
t81 = t50 * qJD(3);
t80 = t53 * qJD(3);
t78 = pkin(7) * t98;
t77 = -0.2e1 * pkin(2) * qJD(3);
t76 = pkin(7) * t81;
t42 = pkin(7) * t80;
t75 = t46 * t88;
t73 = t47 * t80;
t72 = t50 * t80;
t69 = t45 * t47 * t100;
t68 = 0.2e1 * t72;
t43 = t45 ^ 2;
t67 = 0.2e1 * (t47 ^ 2 + t43) * qJD(4);
t66 = (-pkin(7) * t47 + pkin(8)) * t50;
t65 = pkin(4) * t45 - pkin(8) * t47;
t28 = -t48 * t53 + t50 * t97;
t17 = -t28 * qJD(3) + t53 * t74;
t5 = t17 * t45 - t47 * t75;
t6 = t17 * t47 + t45 * t75;
t64 = t5 * t45 + t6 * t47;
t13 = t45 * t76 + t95;
t25 = t45 * t27;
t14 = -t47 * t76 + t25;
t62 = -t13 * t45 + t14 * t47;
t16 = t29 * t47 - t45 * t96;
t61 = t16 * t52 + t28 * t49;
t60 = -t16 * t49 + t28 * t52;
t21 = -t53 * pkin(8) + t24;
t26 = (pkin(7) + t65) * t50;
t59 = t52 * t21 + t49 * t26;
t30 = t47 * t91 + t90;
t34 = -t47 * pkin(4) - t45 * pkin(8) - pkin(3);
t58 = t49 * t34 + t52 * t89;
t57 = t65 * t53;
t56 = -t50 * t84 + t52 * t80;
t55 = t49 * t80 + t50 * t83;
t37 = t45 * t80;
t31 = -t49 * t53 + t52 * t93;
t23 = -t78 + t94;
t15 = t29 * t45 + t47 * t96;
t12 = -t58 * qJD(5) - t47 * t87;
t11 = t102 * t47 - t34 * t83;
t10 = -t71 * t81 - t95;
t9 = t55 * t47 - t52 * t81 - t53 * t84;
t8 = -t30 * qJD(5) + (t47 * t90 + t91) * qJD(3);
t4 = -t49 * t25 + t52 * t42 - t59 * qJD(5) + (-t49 * t66 + t52 * t57) * qJD(3);
t3 = t21 * t84 - t49 * (qJD(3) * t57 + t42) - t26 * t83 - t52 * (qJD(3) * t66 + t25);
t2 = t60 * qJD(5) + t18 * t49 + t6 * t52;
t1 = -t61 * qJD(5) + t18 * t52 - t6 * t49;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t15 * t5 + 0.2e1 * t16 * t6 + 0.2e1 * t28 * t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t75, -t74, 0, 0, 0, 0, 0, (-t53 * t88 - t54 * t81) * t46, (t50 * t88 - t54 * t80) * t46, t45 * t99 + t5 * t53 + (-t15 * t50 + t28 * t98) * qJD(3), t18 * t93 + t6 * t53 + (-t16 * t50 + t28 * t92) * qJD(3), (-t45 * t6 + t47 * t5) * t50 + (t15 * t47 - t16 * t45) * t80, -t15 * t13 + t16 * t14 - t5 * t23 + t6 * t24 + (t28 * t80 + t99) * pkin(7), 0, 0, 0, 0, 0, t15 * t9 + t5 * t30 + (t1 * t50 + t60 * t80) * t45, t15 * t8 + t5 * t31 + (-t2 * t50 - t61 * t80) * t45; 0, 0, 0, 0, t68, 0.2e1 * (-t50 ^ 2 + t53 ^ 2) * qJD(3), 0, 0, 0, t50 * t77, t53 * t77, -0.2e1 * t13 * t53 + 0.2e1 * (t23 + 0.2e1 * t78) * t81, 0.2e1 * t14 * t53 + 0.2e1 * (-t24 + 0.2e1 * t41) * t81, 0.2e1 * (-t13 * t47 - t14 * t45) * t50 + 0.2e1 * (-t23 * t47 - t24 * t45) * t80, 0.2e1 * pkin(7) ^ 2 * t72 + 0.2e1 * t23 * t13 + 0.2e1 * t24 * t14, 0.2e1 * t31 * t8, -0.2e1 * t8 * t30 - 0.2e1 * t31 * t9, (t31 * t80 + t50 * t8) * t101, (-t30 * t80 - t50 * t9) * t101, t43 * t68, 0.2e1 * t10 * t30 + 0.2e1 * t20 * t9 + 0.2e1 * (t4 * t50 + (-t49 * t21 + t52 * t26) * t80) * t45, 0.2e1 * t10 * t31 + 0.2e1 * t20 * t8 + 0.2e1 * (t3 * t50 - t59 * t80) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t17, -t18 * t47, t18 * t45, t64, -t18 * pkin(3) + (t15 * t45 + t16 * t47) * qJD(4) + t64 * qJ(4), 0, 0, 0, 0, 0, -t1 * t47 + (t15 * t83 + t49 * t5) * t45, t2 * t47 + (-t15 * t84 + t5 * t52) * t45; 0, 0, 0, 0, 0, 0, t80, -t81, 0, -t42, t76, t45 * t85 + (t63 * t45 - t41) * qJD(3), t47 * t85 + (t63 * t47 + t78) * qJD(3), t62, -pkin(3) * t42 + (-t23 * t45 + t24 * t47) * qJD(4) + t62 * qJ(4), (-t31 * t84 + t52 * t8) * t45, (-t49 * t8 - t52 * t9 + (t30 * t49 - t31 * t52) * qJD(5)) * t45, t56 * t43 - t8 * t47, -t55 * t43 + t9 * t47, -t45 * t73, -t4 * t47 + (t12 * t50 + (t52 * t34 - t49 * t89) * t80 + qJD(4) * t30 + qJ(4) * t9 + t10 * t49 + t52 * t82) * t45, -t3 * t47 + (qJ(4) * t8 + qJD(4) * t31 + t10 * t52 + t11 * t50 - t49 * t82 - t58 * t80) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, qJ(4) * t67, -0.2e1 * t43 * t49 * t83, (t49 ^ 2 - t52 ^ 2) * t43 * t100, t49 * t69, t52 * t69, 0, -0.2e1 * t12 * t47 + 0.2e1 * (t52 * t79 + t87) * t43, -0.2e1 * t102 * t43 - 0.2e1 * t11 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t73, 0, t42, 0, 0, 0, 0, 0, t56 * t45, -t55 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t84, t47 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, t37, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45 * t84, -t45 * t83, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
