% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRPR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:44
% EndTime: 2019-03-09 01:53:48
% DurationCPUTime: 2.00s
% Computational Cost: add. (2973->194), mult. (6193->356), div. (0->0), fcn. (6220->8), ass. (0->114)
t143 = cos(qJ(6));
t109 = qJD(6) * t143;
t74 = sin(qJ(6));
t120 = qJD(6) * t74;
t69 = sin(pkin(10));
t71 = cos(pkin(10));
t148 = t71 * t109 - t69 * t120;
t142 = sin(qJ(4));
t70 = sin(pkin(9));
t72 = cos(pkin(9));
t75 = cos(qJ(4));
t51 = t142 * t72 + t75 * t70;
t46 = t51 * qJD(4);
t127 = t74 * t71;
t52 = t143 * t69 + t127;
t53 = -t142 * t70 + t75 * t72;
t15 = t148 * t53 - t46 * t52;
t64 = -t71 * pkin(5) - pkin(4);
t151 = 0.2e1 * t64;
t107 = t142 * qJD(3);
t110 = qJD(4) * t142;
t73 = -pkin(1) - qJ(3);
t144 = -pkin(7) + t73;
t114 = t75 * t144;
t119 = t75 * qJD(3);
t54 = t144 * t70;
t149 = (qJD(4) * t114 - t107) * t72 - t54 * t110 - t70 * t119;
t45 = t52 * qJD(6);
t92 = t143 * t71 - t74 * t69;
t13 = t53 * t45 + t92 * t46;
t121 = qJD(4) * t75;
t47 = -t70 * t110 + t72 * t121;
t19 = -t45 * t51 + t92 * t47;
t95 = t46 * t51 - t53 * t47;
t150 = 0.2e1 * t95;
t28 = t92 * t51;
t55 = (t70 ^ 2 + t72 ^ 2) * qJD(3);
t147 = t71 * t149;
t146 = 2 * qJD(2);
t145 = t46 * pkin(4);
t103 = t144 * t142;
t21 = t54 * t121 - t70 * t107 + (qJD(4) * t103 + t119) * t72;
t31 = -t72 * t114 + t142 * t54;
t141 = t31 * t21;
t138 = t46 * t69;
t137 = t46 * t71;
t136 = t92 * t45;
t135 = t51 * t47;
t134 = t52 * t148;
t133 = t52 * t47;
t36 = t53 * t46;
t32 = t72 * t103 + t75 * t54;
t131 = t69 * t32;
t130 = t69 * t47;
t129 = t69 * t53;
t126 = pkin(8) + qJ(5);
t27 = t52 * t53;
t125 = -t148 * t27 - t52 * t15;
t61 = t70 * pkin(3) + qJ(2);
t93 = t51 * pkin(4) + t61;
t88 = -t53 * qJ(5) + t93;
t17 = t71 * t32 + t69 * t88;
t65 = t69 ^ 2;
t67 = t71 ^ 2;
t123 = t65 + t67;
t118 = qJ(2) * qJD(2);
t33 = 0.2e1 * t135;
t117 = -0.2e1 * t36;
t116 = t69 * t137;
t112 = t126 * t69;
t111 = t123 * t47;
t108 = t143 * qJD(5);
t106 = t123 * qJD(5);
t105 = 0.2e1 * t106;
t79 = t149 * t69;
t91 = t47 * pkin(4) - t53 * qJD(5) + qJD(2);
t86 = t46 * qJ(5) + t91;
t8 = t71 * t86 - t79;
t9 = t69 * t86 + t147;
t102 = t9 * t69 + t8 * t71;
t101 = -t8 * t69 + t9 * t71;
t29 = t92 * t53;
t100 = t13 * t92 + t29 * t45;
t16 = t71 * t88 - t131;
t99 = -t16 * t69 + t17 * t71;
t98 = t21 * t53 - t31 * t46;
t97 = -t148 * t51 - t133;
t96 = t36 - t135;
t94 = t143 * t112;
t90 = 0.2e1 * t96;
t89 = -qJ(5) * t47 - qJD(5) * t51 + t145;
t84 = t126 * t46 + t91;
t82 = t51 * pkin(5) - t131 + (-t126 * t53 + t93) * t71;
t81 = t74 * t82;
t80 = t143 * t82;
t78 = t149 * t51 + t32 * t47 - t98;
t77 = t84 * t69 + t147;
t76 = t47 * pkin(5) + t84 * t71 - t79;
t56 = t126 * t71;
t42 = t71 * t47;
t38 = -t74 * t112 + t143 * t56;
t37 = -t74 * t56 - t94;
t26 = t52 * t51;
t24 = -t56 * t109 - qJD(5) * t127 + (t126 * t120 - t108) * t69;
t23 = qJD(6) * t94 - t71 * t108 + (t69 * qJD(5) + qJD(6) * t56) * t74;
t22 = pkin(5) * t129 + t31;
t18 = -pkin(5) * t138 + t21;
t14 = -qJD(6) * t28 - t133;
t11 = -pkin(8) * t129 + t17;
t4 = t143 * t11 + t81;
t3 = -t74 * t11 + t80;
t2 = -qJD(6) * t81 - t11 * t109 + t143 * t76 - t74 * t77;
t1 = -qJD(6) * t80 + t11 * t120 - t143 * t77 - t74 * t76;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, 0.2e1 * t118, 0, 0, 0, 0, 0, 0, t70 * t146, t72 * t146, 0.2e1 * t55, -0.2e1 * t73 * t55 + 0.2e1 * t118, t117, t150, 0, t33, 0, 0, 0.2e1 * qJD(2) * t51 + 0.2e1 * t61 * t47, 0.2e1 * qJD(2) * t53 - 0.2e1 * t61 * t46, -0.2e1 * t78, 0.2e1 * t61 * qJD(2) + 0.2e1 * t149 * t32 + 0.2e1 * t141, t67 * t117, 0.4e1 * t53 * t116, -0.2e1 * t95 * t71, t65 * t117, t69 * t150, t33, 0.2e1 * t16 * t47 + 0.2e1 * t8 * t51 + 0.2e1 * t98 * t69, -0.2e1 * t17 * t47 - 0.2e1 * t9 * t51 + 0.2e1 * t98 * t71, -0.2e1 * t102 * t53 + 0.2e1 * (t16 * t71 + t17 * t69) * t46, 0.2e1 * t16 * t8 + 0.2e1 * t17 * t9 + 0.2e1 * t141, -0.2e1 * t29 * t13, 0.2e1 * t13 * t27 - 0.2e1 * t29 * t15, -0.2e1 * t13 * t51 + 0.2e1 * t29 * t47, 0.2e1 * t27 * t15, -0.2e1 * t15 * t51 - 0.2e1 * t27 * t47, t33, 0.2e1 * t22 * t15 + 0.2e1 * t18 * t27 + 0.2e1 * t2 * t51 + 0.2e1 * t3 * t47, 0.2e1 * t1 * t51 - 0.2e1 * t22 * t13 + 0.2e1 * t18 * t29 - 0.2e1 * t4 * t47, 0.2e1 * t1 * t27 + 0.2e1 * t3 * t13 - 0.2e1 * t4 * t15 - 0.2e1 * t2 * t29, -0.2e1 * t4 * t1 + 0.2e1 * t22 * t18 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, 0, 0, 0, 0, 0, 0, 0, 0, t90, t78, 0, 0, 0, 0, 0, 0, t69 * t90, t71 * t90, 0, t101 * t51 + t99 * t47 - t98, 0, 0, 0, 0, 0, 0, t14 * t51 - t53 * t15 - t26 * t47 + t46 * t27, t53 * t13 - t19 * t51 - t28 * t47 + t46 * t29, -t26 * t13 - t14 * t29 - t28 * t15 - t19 * t27, -t1 * t28 + t3 * t14 - t18 * t53 + t19 * t4 - t2 * t26 + t22 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t51 * t111 - 0.2e1 * t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t26 * t14 + 0.2e1 * t19 * t28 - 0.2e1 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, t47, -t46, 0, qJD(2), 0, 0, 0, 0, 0, 0, t42, -t130, t123 * t46, t102, 0, 0, 0, 0, 0, 0, t19, t97, t100 + t125, -t1 * t52 + t148 * t4 + t2 * t92 - t3 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t92 + t148 * t28 + t19 * t52 + t26 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t134 - 0.2e1 * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, -t47, 0, -t21, -t149, 0, 0, -t116 (t65 - t67) * t46, t130, t116, t42, 0, -t21 * t71 + t89 * t69, t21 * t69 + t89 * t71, t101, -t21 * pkin(4) + t101 * qJ(5) + t99 * qJD(5), -t13 * t52 + t148 * t29, -t100 + t125, -t97, -t15 * t92 + t27 * t45, t19, 0, t64 * t15 - t18 * t92 + t22 * t45 + t24 * t51 + t37 * t47, -t64 * t13 + t148 * t22 + t18 * t52 + t23 * t51 - t38 * t47, -t1 * t92 + t37 * t13 - t148 * t3 - t38 * t15 - t2 * t52 + t23 * t27 - t24 * t29 - t4 * t45, -t1 * t38 + t18 * t64 + t2 * t37 - t4 * t23 + t3 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t47, 0, 0, 0, 0, 0, 0, 0, 0, -t137, t138, t111, qJ(5) * t111 + t51 * t106 - t145, 0, 0, 0, 0, 0, 0, -t13, -t15, -t14 * t52 + t148 * t26 + t19 * t92 - t28 * t45, t14 * t37 + t19 * t38 - t28 * t23 - t26 * t24 + t46 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148 * t38 - t52 * t23 + t24 * t92 - t45 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, qJ(5) * t105, 0.2e1 * t134, 0.2e1 * t148 * t92 - 0.2e1 * t52 * t45, 0, -0.2e1 * t136, 0, 0, t45 * t151, t148 * t151, -0.2e1 * t148 * t37 - 0.2e1 * t23 * t92 - 0.2e1 * t24 * t52 - 0.2e1 * t38 * t45, -0.2e1 * t38 * t23 + 0.2e1 * t37 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, -t137, 0, t21, 0, 0, 0, 0, 0, 0, t15, -t13, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t148, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, -t15, t47, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t148, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, 0, -t45, 0, t24, t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
