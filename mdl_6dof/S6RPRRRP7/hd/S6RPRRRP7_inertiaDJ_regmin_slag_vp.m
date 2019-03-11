% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:21:10
% EndTime: 2019-03-09 06:21:15
% DurationCPUTime: 2.06s
% Computational Cost: add. (4317->257), mult. (9672->417), div. (0->0), fcn. (9559->8), ass. (0->129)
t105 = sin(qJ(4));
t107 = cos(qJ(4));
t148 = qJD(4) * t107;
t102 = sin(pkin(10));
t103 = cos(pkin(10));
t106 = sin(qJ(3));
t151 = t106 * t103;
t166 = cos(qJ(3));
t78 = t166 * t102 + t151;
t140 = t78 * t148;
t111 = -t106 * t102 + t166 * t103;
t65 = t111 * qJD(3);
t173 = t105 * t65 + t140;
t165 = cos(qJ(5));
t131 = t165 * qJD(5);
t172 = t165 * qJD(4) + t131;
t137 = t165 * t107;
t104 = sin(qJ(5));
t154 = t104 * t105;
t79 = -t137 + t154;
t101 = t107 ^ 2;
t150 = t105 ^ 2 - t101;
t128 = t150 * qJD(4);
t171 = qJD(4) + qJD(5);
t167 = t111 * pkin(4);
t168 = pkin(9) * t78;
t93 = -t103 * pkin(2) - pkin(1);
t46 = -pkin(3) * t111 - t78 * pkin(8) + t93;
t43 = t107 * t46;
t162 = pkin(7) + qJ(2);
t81 = t162 * t102;
t156 = t106 * t81;
t82 = t162 * t103;
t53 = t166 * t82 - t156;
t21 = -t105 * t53 - t107 * t168 - t167 + t43;
t157 = t105 * t78;
t48 = t107 * t53;
t159 = t105 * t46 + t48;
t23 = -pkin(9) * t157 + t159;
t22 = t165 * t23;
t114 = t104 * t21 + t22;
t149 = qJD(4) * t105;
t133 = qJD(3) * t166;
t134 = qJD(2) * t166;
t32 = t81 * t133 - t103 * t134 + (qJD(2) * t102 + qJD(3) * t82) * t106;
t66 = t78 * qJD(3);
t45 = t66 * pkin(3) - t65 * pkin(8);
t13 = -t105 * t45 + t107 * t32 - t46 * t148 + t53 * t149;
t12 = -pkin(9) * t173 - t13;
t130 = t105 * t32 + t107 * t45;
t155 = t107 * t65;
t8 = -pkin(9) * t155 + t66 * pkin(4) + (-t48 + (-t46 + t168) * t105) * qJD(4) + t130;
t135 = t104 * t12 - t165 * t8;
t4 = -t114 * qJD(5) - t135;
t170 = 0.2e1 * t93;
t108 = 2 * qJD(6);
t169 = -pkin(9) - pkin(8);
t63 = t66 * pkin(5);
t164 = t78 * t65;
t50 = -t107 * t172 + t154 * t171;
t138 = t165 * t105;
t153 = t104 * t107;
t80 = t138 + t153;
t163 = t80 * t50;
t141 = t78 * t149;
t142 = t78 * t154;
t17 = t65 * t138 - t104 * t141 - qJD(5) * t142 + (t104 * t65 + t172 * t78) * t107;
t40 = t80 * t78;
t161 = -t80 * t17 + t50 * t40;
t152 = t105 * t107;
t147 = qJD(5) * t104;
t146 = -0.2e1 * pkin(3) * qJD(4);
t145 = t165 * pkin(4);
t144 = pkin(4) * t149;
t143 = pkin(4) * t147;
t96 = -t107 * pkin(4) - pkin(3);
t139 = t169 * qJD(4);
t136 = t105 * t148;
t129 = -0.4e1 * t78 * t152;
t33 = qJD(2) * t151 - qJD(3) * t156 + t102 * t134 + t82 * t133;
t127 = 0.2e1 * (t102 ^ 2 + t103 ^ 2) * qJD(2);
t126 = t169 * t165;
t125 = -pkin(3) * t65 - pkin(8) * t66;
t124 = pkin(3) * t78 - pkin(8) * t111;
t52 = t106 * t82 + t166 * t81;
t51 = t171 * t80;
t16 = t51 * t78 + t65 * t79;
t41 = t78 * t137 - t142;
t123 = -t79 * t16 + t51 * t41;
t116 = qJD(4) * t126;
t117 = t105 * t126;
t83 = t169 * t107;
t29 = -qJD(5) * t117 - t105 * t116 - t139 * t153 - t83 * t147;
t55 = t169 * t154 - t165 * t83;
t122 = t111 * t29 + t55 * t66;
t30 = -t83 * t131 - t107 * t116 + (qJD(5) * t169 + t139) * t154;
t54 = -t104 * t83 - t117;
t121 = t111 * t30 - t54 * t66;
t120 = t33 * t78 + t52 * t65;
t119 = -t111 * t50 - t80 * t66;
t27 = t111 * t51 - t79 * t66;
t118 = -t111 * t65 + t78 * t66;
t24 = t173 * pkin(4) + t33;
t34 = pkin(4) * t157 + t52;
t115 = -t104 * t23 + t165 * t21;
t112 = t105 * t66 - t111 * t148;
t3 = -t104 * t8 - t165 * t12 - t21 * t131 + t23 * t147;
t110 = -t51 * pkin(5) - t50 * qJ(6) + t80 * qJD(6);
t62 = t66 * qJ(6);
t67 = t111 * qJD(6);
t1 = -t3 + t62 - t67;
t109 = (-t22 + (-t21 + t167) * t104) * qJD(5) - t135;
t97 = pkin(4) * t131;
t95 = -t145 - pkin(5);
t92 = t104 * pkin(4) + qJ(6);
t90 = -0.2e1 * t143;
t84 = t97 + qJD(6);
t74 = t78 ^ 2;
t49 = -0.2e1 * t111 * t66;
t47 = t79 * pkin(5) - t80 * qJ(6) + t96;
t44 = t107 * t66 + t111 * t149;
t25 = -t110 + t144;
t18 = t40 * pkin(5) - t41 * qJ(6) + t34;
t14 = -t159 * qJD(4) + t130;
t10 = pkin(5) * t111 - t115;
t9 = -qJ(6) * t111 + t114;
t5 = t17 * pkin(5) + t16 * qJ(6) - t41 * qJD(6) + t24;
t2 = -t4 - t63;
t6 = [0, 0, 0, 0, 0, t127, qJ(2) * t127, 0.2e1 * t164, -0.2e1 * t118, 0, 0, 0, t66 * t170, t65 * t170, 0.2e1 * t101 * t164 - 0.2e1 * t74 * t136, 0.2e1 * t128 * t74 + t65 * t129, 0.2e1 * t107 * t118 + 0.2e1 * t111 * t141, -0.2e1 * t105 * t118 + 0.2e1 * t111 * t140, t49, 0.2e1 * t52 * t140 - 0.2e1 * t14 * t111 + 0.2e1 * t43 * t66 + 0.2e1 * (-t53 * t66 + t120) * t105, 0.2e1 * t107 * t120 - 0.2e1 * t111 * t13 - 0.2e1 * t141 * t52 - 0.2e1 * t159 * t66, -0.2e1 * t41 * t16, 0.2e1 * t16 * t40 - 0.2e1 * t41 * t17, 0.2e1 * t111 * t16 + 0.2e1 * t41 * t66, 0.2e1 * t111 * t17 - 0.2e1 * t40 * t66, t49, -0.2e1 * t111 * t4 + 0.2e1 * t115 * t66 + 0.2e1 * t34 * t17 + 0.2e1 * t24 * t40, -0.2e1 * t111 * t3 - 0.2e1 * t114 * t66 - 0.2e1 * t34 * t16 + 0.2e1 * t24 * t41, -0.2e1 * t10 * t66 + 0.2e1 * t111 * t2 + 0.2e1 * t18 * t17 + 0.2e1 * t5 * t40, -0.2e1 * t1 * t40 - 0.2e1 * t10 * t16 - 0.2e1 * t9 * t17 + 0.2e1 * t2 * t41, -0.2e1 * t1 * t111 + 0.2e1 * t18 * t16 - 0.2e1 * t5 * t41 + 0.2e1 * t9 * t66, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t18 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t65, 0, 0, 0, 0, 0, t44, -t112, 0, 0, 0, 0, 0, t27, t119, t27, t123 + t161, -t119, t1 * t80 + t10 * t51 + t2 * t79 - t9 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t79 * t51 - 0.2e1 * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t66, 0, -t33, t32, -t78 * t128 + t65 * t152, qJD(4) * t129 - t150 * t65, t112, t44, 0, -t33 * t107 + t125 * t105 + (t105 * t52 - t107 * t124) * qJD(4), t33 * t105 + t125 * t107 + (t105 * t124 + t107 * t52) * qJD(4), -t16 * t80 - t41 * t50, -t123 + t161, -t119, t27, 0, t144 * t40 + t96 * t17 + t24 * t79 + t34 * t51 + t121, t144 * t41 - t96 * t16 + t24 * t80 - t34 * t50 - t122, t47 * t17 + t18 * t51 + t25 * t40 + t5 * t79 + t121, -t1 * t79 - t10 * t50 - t54 * t16 - t55 * t17 + t2 * t80 + t29 * t40 + t30 * t41 - t9 * t51, t47 * t16 + t18 * t50 - t25 * t41 - t5 * t80 + t122, t1 * t55 + t10 * t30 + t18 * t25 + t2 * t54 - t9 * t29 + t5 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t29 + t79 * t30 - t50 * t55 + t51 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t136, -0.2e1 * t128, 0, 0, 0, t105 * t146, t107 * t146, -0.2e1 * t163, 0.2e1 * t79 * t50 - 0.2e1 * t51 * t80, 0, 0, 0, 0.2e1 * t144 * t79 + 0.2e1 * t96 * t51, 0.2e1 * t144 * t80 - 0.2e1 * t96 * t50, 0.2e1 * t25 * t79 + 0.2e1 * t47 * t51, 0.2e1 * t29 * t79 + 0.2e1 * t30 * t80 - 0.2e1 * t54 * t50 - 0.2e1 * t55 * t51, -0.2e1 * t25 * t80 + 0.2e1 * t47 * t50, 0.2e1 * t47 * t25 - 0.2e1 * t55 * t29 + 0.2e1 * t54 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141 + t155, -t173, t66, t14, t13, 0, 0, -t16, -t17, t66, t145 * t66 + t109 (-t104 * t66 + t111 * t131) * pkin(4) + t3, -t95 * t66 + t109 + t63, t143 * t41 - t95 * t16 - t92 * t17 - t84 * t40, -t111 * t84 + t92 * t66 + t1, t1 * t92 + t10 * t143 + t2 * t95 + t9 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, -t148, 0, 0, 0, 0, 0, -t51, t50, -t51, 0, -t50, t143 * t79 - t50 * t92 + t51 * t95 + t80 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t149, 0, -pkin(8) * t148, pkin(8) * t149, 0, 0, -t50, -t51, 0, -t30, t29, -t30, t143 * t80 - t95 * t50 - t92 * t51 - t84 * t79, -t29, t143 * t54 - t29 * t92 + t30 * t95 + t55 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -0.2e1 * t97, t90, 0, 0.2e1 * t84, 0.2e1 * t143 * t95 + 0.2e1 * t92 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, t66, t4, t3, t4 + 0.2e1 * t63, pkin(5) * t16 - t17 * qJ(6) - t40 * qJD(6), -t3 + 0.2e1 * t62 - 0.2e1 * t67, -t2 * pkin(5) + t1 * qJ(6) + t9 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50, -t51, 0, -t50, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t51, 0, -t30, t29, -t30, pkin(5) * t50 - t51 * qJ(6) - t79 * qJD(6), -t29, -t30 * pkin(5) - t29 * qJ(6) + t55 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, -t97, -t143, 0, t108 + t97, -pkin(5) * t143 + t84 * qJ(6) + t92 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, qJ(6) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t16, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
