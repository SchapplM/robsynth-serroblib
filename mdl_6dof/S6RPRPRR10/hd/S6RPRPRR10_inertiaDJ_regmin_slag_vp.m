% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:09:20
% EndTime: 2019-03-09 04:09:23
% DurationCPUTime: 1.46s
% Computational Cost: add. (1903->207), mult. (4535->379), div. (0->0), fcn. (4303->8), ass. (0->121)
t92 = cos(pkin(10));
t97 = cos(qJ(5));
t140 = t97 * t92;
t91 = sin(pkin(10));
t94 = sin(qJ(5));
t71 = t91 * t94 - t140;
t95 = sin(qJ(3));
t156 = t71 * t95;
t72 = t91 * t97 + t92 * t94;
t98 = cos(qJ(3));
t155 = t72 * t98;
t99 = -pkin(1) - pkin(7);
t118 = -t91 * t99 + pkin(4);
t149 = pkin(8) * t92;
t148 = t95 * pkin(3);
t111 = -qJ(4) * t98 + t148;
t73 = qJ(2) + t111;
t66 = t92 * t73;
t40 = t118 * t95 - t98 * t149 + t66;
t144 = t91 * t98;
t142 = t95 * t99;
t81 = t92 * t142;
t49 = t91 * t73 + t81;
t47 = -pkin(8) * t144 + t49;
t138 = t94 * t40 + t97 * t47;
t139 = pkin(8) + qJ(4);
t78 = t139 * t91;
t79 = t139 * t92;
t136 = -t94 * t78 + t97 * t79;
t154 = t98 * t140 - t94 * t144;
t132 = qJD(5) * t97;
t24 = (qJD(4) * t91 + qJD(5) * t79) * t94 - qJD(4) * t140 + t78 * t132;
t87 = -pkin(4) * t92 - pkin(3);
t153 = 0.2e1 * t87;
t152 = 2 * qJD(2);
t63 = t72 * qJD(5);
t151 = pkin(5) * t63;
t150 = pkin(5) * t95;
t131 = t95 * qJD(3);
t122 = t92 * t131;
t123 = t91 * t131;
t35 = t154 * qJD(5) - t94 * t122 - t97 * t123;
t146 = t47 * t94;
t59 = -t98 * qJD(4) + qJD(2) + (pkin(3) * t98 + qJ(4) * t95) * qJD(3);
t51 = t92 * t59;
t23 = t51 + (t118 * t98 + t95 * t149) * qJD(3);
t88 = t98 * qJD(3);
t120 = t99 * t88;
t46 = t92 * t120 + t91 * t59;
t36 = pkin(8) * t123 + t46;
t8 = qJD(5) * t146 - t40 * t132 - t94 * t23 - t97 * t36;
t7 = -pkin(9) * t35 - t8;
t96 = cos(qJ(6));
t147 = t96 * t7;
t145 = t63 * t95;
t15 = -pkin(9) * t155 + t138;
t93 = sin(qJ(6));
t143 = t93 * t15;
t141 = t96 * t15;
t135 = t91 ^ 2 + t92 ^ 2;
t134 = pkin(5) * qJD(6);
t133 = qJD(4) * t95;
t130 = qJ(2) * qJD(3);
t128 = t91 * t142;
t126 = pkin(5) * t88;
t125 = t93 * t134;
t124 = t96 * t134;
t85 = t99 * t131;
t121 = t95 * t88;
t33 = qJD(3) * t156 - qJD(5) * t155;
t9 = -t138 * qJD(5) + t97 * t23 - t94 * t36;
t6 = -t33 * pkin(9) + t126 + t9;
t119 = t96 * t6 - t93 * t7;
t115 = t97 * t40 - t146;
t58 = t71 * t98;
t14 = pkin(9) * t58 + t115 + t150;
t117 = -t14 - t150;
t116 = t135 * t98;
t114 = -t97 * t78 - t79 * t94;
t67 = pkin(4) * t144 - t98 * t99;
t113 = t135 * qJD(4);
t82 = 0.2e1 * t121;
t112 = 0.2e1 * t113;
t110 = t14 * t96 - t143;
t109 = t14 * t93 + t141;
t30 = -pkin(9) * t72 + t114;
t31 = -pkin(9) * t71 + t136;
t108 = t30 * t96 - t31 * t93;
t107 = t30 * t93 + t31 * t96;
t45 = -t91 * t120 + t51;
t106 = -t45 * t91 + t46 * t92;
t48 = t66 - t128;
t105 = -t48 * t91 + t49 * t92;
t55 = t72 * t95;
t104 = t156 * t93 - t55 * t96;
t103 = -t156 * t96 - t55 * t93;
t28 = t155 * t96 - t58 * t93;
t29 = -t155 * t93 - t58 * t96;
t42 = t96 * t71 + t72 * t93;
t43 = -t71 * t93 + t72 * t96;
t61 = -pkin(4) * t123 + t85;
t25 = -t72 * qJD(4) - t136 * qJD(5);
t62 = t71 * qJD(5);
t53 = pkin(5) * t71 + t87;
t44 = pkin(5) * t155 + t67;
t34 = -qJD(3) * t155 + qJD(5) * t156;
t32 = -t154 * qJD(3) + t145;
t20 = pkin(5) * t35 + t61;
t19 = t62 * pkin(9) + t25;
t18 = -t63 * pkin(9) - t24;
t17 = t43 * qJD(6) - t93 * t62 + t96 * t63;
t16 = -t42 * qJD(6) - t96 * t62 - t93 * t63;
t13 = t29 * qJD(6) + t93 * t33 + t96 * t35;
t12 = -t103 * qJD(6) + t93 * t32 + t96 * t34;
t11 = -t28 * qJD(6) + t96 * t33 - t93 * t35;
t10 = -t104 * qJD(6) + t96 * t32 - t93 * t34;
t4 = -t107 * qJD(6) - t93 * t18 + t96 * t19;
t3 = -t108 * qJD(6) - t96 * t18 - t93 * t19;
t2 = -t109 * qJD(6) + t119;
t1 = -t110 * qJD(6) - t93 * t6 - t147;
t5 = [0, 0, 0, 0, t152, qJ(2) * t152, -0.2e1 * t121, 0.2e1 * (t95 ^ 2 - t98 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t95 + 0.2e1 * t98 * t130, 0.2e1 * qJD(2) * t98 - 0.2e1 * t95 * t130, 0.2e1 * t45 * t95 + 0.2e1 * (t48 + 0.2e1 * t128) * t88, -0.2e1 * t46 * t95 + 0.2e1 * (-t49 + 0.2e1 * t81) * t88, 0.2e1 * (-t45 * t92 - t46 * t91) * t98 + 0.2e1 * (t48 * t92 + t49 * t91) * t131, -0.2e1 * t99 ^ 2 * t121 + 0.2e1 * t48 * t45 + 0.2e1 * t49 * t46, -0.2e1 * t58 * t33, -0.2e1 * t155 * t33 + 0.2e1 * t35 * t58, 0.2e1 * t33 * t95 - 0.2e1 * t58 * t88, -0.2e1 * t155 * t88 - 0.2e1 * t35 * t95, t82, 0.2e1 * t115 * t88 + 0.2e1 * t155 * t61 + 0.2e1 * t67 * t35 + 0.2e1 * t9 * t95, -0.2e1 * t138 * t88 + 0.2e1 * t67 * t33 - 0.2e1 * t61 * t58 + 0.2e1 * t8 * t95, 0.2e1 * t29 * t11, -0.2e1 * t11 * t28 - 0.2e1 * t13 * t29, 0.2e1 * t11 * t95 + 0.2e1 * t29 * t88, -0.2e1 * t13 * t95 - 0.2e1 * t28 * t88, t82, 0.2e1 * t110 * t88 + 0.2e1 * t44 * t13 + 0.2e1 * t2 * t95 + 0.2e1 * t20 * t28, 0.2e1 * t1 * t95 - 0.2e1 * t109 * t88 + 0.2e1 * t44 * t11 + 0.2e1 * t20 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106 * t95 + (t105 - 0.2e1 * t142) * t88, 0, 0, 0, 0, 0, t34 * t95 - t98 * t35 + (t155 * t95 - t55 * t98) * qJD(3), t32 * t95 - t98 * t33 + (t156 * t98 - t58 * t95) * qJD(3), 0, 0, 0, 0, 0, t12 * t95 - t98 * t13 + (t104 * t98 + t95 * t28) * qJD(3), t10 * t95 - t98 * t11 + (-t103 * t98 + t95 * t29) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t135) * t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t88, 0, -t85, -t120, -t91 * t133 + (t111 * t91 - t81) * qJD(3), -t92 * t133 + (t111 * t92 + t128) * qJD(3), t106, -pkin(3) * t85 + t106 * qJ(4) + t105 * qJD(4), t33 * t72 + t58 * t62, t155 * t62 - t33 * t71 - t35 * t72 + t58 * t63, -t62 * t95 + t72 * t88, -t71 * t88 - t145, 0, t114 * t88 + t25 * t95 + t87 * t35 + t61 * t71 + t67 * t63, -t136 * t88 + t24 * t95 + t87 * t33 + t61 * t72 - t67 * t62, t11 * t43 + t16 * t29, -t11 * t42 - t13 * t43 - t16 * t28 - t17 * t29, t16 * t95 + t43 * t88, -t17 * t95 - t42 * t88, 0, t108 * t88 + t53 * t13 + t28 * t151 + t44 * t17 + t20 * t42 + t4 * t95, -t107 * t88 + t53 * t11 + t29 * t151 + t44 * t16 + t20 * t43 + t3 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t88, -t122, t123, qJD(3) * t116, t95 * t113 + (qJ(4) * t116 - t148) * qJD(3), 0, 0, 0, 0, 0, t71 * t131 - t63 * t98, t72 * t131 + t62 * t98, 0, 0, 0, 0, 0, t131 * t42 - t17 * t98, t131 * t43 - t16 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, qJ(4) * t112, -0.2e1 * t72 * t62, 0.2e1 * t62 * t71 - 0.2e1 * t63 * t72, 0, 0, 0, t63 * t153, -t62 * t153, 0.2e1 * t43 * t16, -0.2e1 * t16 * t42 - 0.2e1 * t17 * t43, 0, 0, 0, 0.2e1 * t42 * t151 + 0.2e1 * t17 * t53, 0.2e1 * t43 * t151 + 0.2e1 * t16 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t122, 0, t85, 0, 0, 0, 0, 0, t35, t33, 0, 0, 0, 0, 0, t13, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t62, 0, 0, 0, 0, 0, t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t35, t88, t9, t8, 0, 0, t11, -t13, t88, t96 * t126 + (t117 * t93 - t141) * qJD(6) + t119, -t147 + (-t6 - t126) * t93 + (t117 * t96 + t143) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t32, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t63, 0, t25, t24, 0, 0, t16, -t17, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t125, -0.2e1 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t13, t88, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
