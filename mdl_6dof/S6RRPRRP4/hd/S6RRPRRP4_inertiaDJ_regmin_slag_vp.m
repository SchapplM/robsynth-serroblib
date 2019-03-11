% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:39
% EndTime: 2019-03-09 11:55:47
% DurationCPUTime: 2.22s
% Computational Cost: add. (4571->264), mult. (10024->440), div. (0->0), fcn. (9676->8), ass. (0->136)
t101 = sin(qJ(4));
t103 = cos(qJ(4));
t151 = qJD(4) * t103;
t102 = sin(qJ(2));
t104 = cos(qJ(2));
t157 = sin(pkin(10));
t99 = cos(pkin(10));
t78 = t102 * t99 + t104 * t157;
t142 = t78 * t151;
t131 = t157 * t102;
t153 = qJD(2) * t104;
t73 = -qJD(2) * t131 + t153 * t99;
t173 = t101 * t73 + t142;
t167 = cos(qJ(5));
t135 = t167 * qJD(5);
t172 = qJD(4) * t167 + t135;
t139 = t167 * t103;
t100 = sin(qJ(5));
t156 = t100 * t101;
t80 = -t139 + t156;
t171 = qJD(4) + qJD(5);
t77 = -t104 * t99 + t131;
t168 = t77 * pkin(4);
t169 = pkin(9) * t78;
t141 = -pkin(2) * t104 - pkin(1);
t49 = pkin(3) * t77 - pkin(8) * t78 + t141;
t45 = t103 * t49;
t165 = -qJ(3) - pkin(7);
t83 = t165 * t102;
t84 = t165 * t104;
t55 = t157 * t83 - t84 * t99;
t21 = -t101 * t55 - t103 * t169 + t168 + t45;
t159 = t101 * t78;
t51 = t103 * t55;
t162 = t101 * t49 + t51;
t23 = -pkin(9) * t159 + t162;
t22 = t167 * t23;
t120 = t100 * t21 + t22;
t152 = qJD(4) * t101;
t133 = qJD(2) * t165;
t70 = t104 * qJD(3) + t102 * t133;
t71 = -t102 * qJD(3) + t104 * t133;
t38 = t157 * t71 + t70 * t99;
t72 = t78 * qJD(2);
t154 = qJD(2) * t102;
t96 = pkin(2) * t154;
t39 = pkin(3) * t72 - pkin(8) * t73 + t96;
t13 = -t101 * t39 - t103 * t38 - t151 * t49 + t152 * t55;
t12 = -pkin(9) * t173 - t13;
t134 = -t101 * t38 + t103 * t39;
t158 = t103 * t73;
t8 = -pkin(9) * t158 + t72 * pkin(4) + (-t51 + (-t49 + t169) * t101) * qJD(4) + t134;
t137 = t100 * t12 - t167 * t8;
t4 = -qJD(5) * t120 - t137;
t170 = 0.2e1 * qJD(4);
t105 = 2 * qJD(6);
t68 = t72 * pkin(5);
t52 = -t103 * t172 + t156 * t171;
t140 = t167 * t101;
t81 = t100 * t103 + t140;
t166 = t81 * t52;
t143 = t78 * t152;
t144 = t78 * t156;
t17 = t73 * t140 - t100 * t143 - qJD(5) * t144 + (t100 * t73 + t172 * t78) * t103;
t41 = t81 * t78;
t164 = -t17 * t81 + t41 * t52;
t98 = t103 ^ 2;
t161 = t101 ^ 2 - t98;
t155 = t101 * t103;
t150 = qJD(5) * t100;
t149 = -0.2e1 * pkin(1) * qJD(2);
t92 = -pkin(2) * t99 - pkin(3);
t148 = t92 * t170;
t147 = t167 * pkin(4);
t146 = pkin(4) * t152;
t145 = pkin(4) * t150;
t138 = t101 * t151;
t132 = t161 * qJD(4);
t130 = -0.4e1 * t78 * t155;
t129 = pkin(2) * t157 + pkin(8);
t53 = t171 * t81;
t16 = t53 * t78 + t73 * t80;
t42 = t139 * t78 - t144;
t128 = -t16 * t80 + t42 * t53;
t122 = pkin(9) + t129;
t112 = t122 * t167;
t107 = qJD(4) * t112;
t108 = t101 * t112;
t116 = t100 * t122;
t111 = qJD(4) * t116;
t75 = t122 * t103;
t26 = qJD(5) * t108 + t101 * t107 + t103 * t111 + t150 * t75;
t113 = t101 * t116;
t48 = t167 * t75 - t113;
t127 = -t26 * t77 + t48 * t72;
t27 = -qJD(5) * t113 - t101 * t111 + t103 * t107 + t135 * t75;
t47 = t100 * t75 + t108;
t126 = -t27 * t77 - t47 * t72;
t37 = t157 * t70 - t71 * t99;
t54 = -t157 * t84 - t83 * t99;
t125 = t37 * t78 + t54 * t73;
t124 = t52 * t77 - t72 * t81;
t29 = -t53 * t77 - t72 * t80;
t123 = t72 * t78 + t73 * t77;
t82 = -pkin(4) * t103 + t92;
t34 = pkin(4) * t159 + t54;
t121 = -t100 * t23 + t167 * t21;
t119 = t101 * t72 + t151 * t77;
t3 = -t100 * t8 - t12 * t167 - t135 * t21 + t150 * t23;
t117 = qJD(4) * t129;
t25 = pkin(4) * t173 + t37;
t115 = -t55 * t72 + t125;
t114 = -pkin(5) * t53 - qJ(6) * t52 + qJD(6) * t81;
t67 = t72 * qJ(6);
t74 = t77 * qJD(6);
t1 = -t3 + t67 + t74;
t110 = -t129 * t72 + t92 * t73;
t109 = t129 * t77 - t92 * t78;
t106 = (-t22 + (-t21 - t168) * t100) * qJD(5) - t137;
t95 = pkin(4) * t135;
t94 = -t147 - pkin(5);
t91 = pkin(4) * t100 + qJ(6);
t89 = -0.2e1 * t145;
t85 = t95 + qJD(6);
t76 = t78 ^ 2;
t50 = 0.2e1 * t77 * t72;
t46 = t103 * t72 - t152 * t77;
t43 = pkin(5) * t80 - qJ(6) * t81 + t82;
t24 = -t114 + t146;
t18 = pkin(5) * t41 - qJ(6) * t42 + t34;
t14 = -qJD(4) * t162 + t134;
t10 = -pkin(5) * t77 - t121;
t9 = qJ(6) * t77 + t120;
t5 = pkin(5) * t17 + qJ(6) * t16 - qJD(6) * t42 + t25;
t2 = -t4 - t68;
t6 = [0, 0, 0, 0.2e1 * t102 * t153, 0.2e1 * (-t102 ^ 2 + t104 ^ 2) * qJD(2), 0, 0, 0, t102 * t149, t104 * t149, -0.2e1 * t38 * t77 + 0.2e1 * t115, 0.2e1 * t141 * t96 + 0.2e1 * t54 * t37 + 0.2e1 * t55 * t38, 0.2e1 * t73 * t78 * t98 - 0.2e1 * t138 * t76, t161 * t170 * t76 + t130 * t73, 0.2e1 * t103 * t123 - 0.2e1 * t143 * t77, -0.2e1 * t101 * t123 - 0.2e1 * t142 * t77, t50, 0.2e1 * t101 * t115 + 0.2e1 * t14 * t77 + 0.2e1 * t142 * t54 + 0.2e1 * t45 * t72, 0.2e1 * t103 * t125 + 0.2e1 * t13 * t77 - 0.2e1 * t143 * t54 - 0.2e1 * t162 * t72, -0.2e1 * t42 * t16, 0.2e1 * t16 * t41 - 0.2e1 * t17 * t42, -0.2e1 * t16 * t77 + 0.2e1 * t42 * t72, -0.2e1 * t17 * t77 - 0.2e1 * t41 * t72, t50, 0.2e1 * t121 * t72 + 0.2e1 * t17 * t34 + 0.2e1 * t25 * t41 + 0.2e1 * t4 * t77, -0.2e1 * t120 * t72 - 0.2e1 * t16 * t34 + 0.2e1 * t25 * t42 + 0.2e1 * t3 * t77, -0.2e1 * t10 * t72 + 0.2e1 * t17 * t18 - 0.2e1 * t2 * t77 + 0.2e1 * t41 * t5, -0.2e1 * t1 * t41 - 0.2e1 * t10 * t16 - 0.2e1 * t17 * t9 + 0.2e1 * t2 * t42, 0.2e1 * t1 * t77 + 0.2e1 * t16 * t18 - 0.2e1 * t42 * t5 + 0.2e1 * t72 * t9, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t18 * t5; 0, 0, 0, 0, 0, t153, -t154, 0, -pkin(7) * t153, pkin(7) * t154 (-t157 * t72 - t73 * t99) * pkin(2) (t157 * t38 - t37 * t99) * pkin(2), -t132 * t78 + t155 * t73, qJD(4) * t130 - t161 * t73, t119, t46, 0, -t37 * t103 + t110 * t101 + (t54 * t101 - t103 * t109) * qJD(4), t37 * t101 + t110 * t103 + (t101 * t109 + t54 * t103) * qJD(4), -t16 * t81 - t42 * t52, -t128 + t164, -t124, t29, 0, t146 * t41 + t17 * t82 + t25 * t80 + t34 * t53 + t126, t146 * t42 - t16 * t82 + t25 * t81 - t34 * t52 - t127, t17 * t43 + t18 * t53 + t24 * t41 + t5 * t80 + t126, -t1 * t80 - t10 * t52 - t16 * t47 - t17 * t48 + t2 * t81 + t26 * t41 + t27 * t42 - t53 * t9, t16 * t43 + t18 * t52 - t24 * t42 - t5 * t81 + t127, t1 * t48 + t10 * t27 + t18 * t24 + t2 * t47 - t26 * t9 + t43 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t138, -0.2e1 * t132, 0, 0, 0, t101 * t148, t103 * t148, -0.2e1 * t166, 0.2e1 * t52 * t80 - 0.2e1 * t53 * t81, 0, 0, 0, 0.2e1 * t146 * t80 + 0.2e1 * t53 * t82, 0.2e1 * t146 * t81 - 0.2e1 * t52 * t82, 0.2e1 * t24 * t80 + 0.2e1 * t43 * t53, 0.2e1 * t26 * t80 + 0.2e1 * t27 * t81 - 0.2e1 * t47 * t52 - 0.2e1 * t48 * t53, -0.2e1 * t24 * t81 + 0.2e1 * t43 * t52, 0.2e1 * t24 * t43 - 0.2e1 * t26 * t48 + 0.2e1 * t27 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, 0, 0, 0, 0, t46, -t119, 0, 0, 0, 0, 0, t29, t124, t29, t128 + t164, -t124, t1 * t81 + t10 * t53 + t2 * t80 - t52 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26 * t81 + t27 * t80 + t47 * t53 - t48 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t53 * t80 - 0.2e1 * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143 + t158, -t173, t72, t14, t13, 0, 0, -t16, -t17, t72, t147 * t72 + t106 (-t100 * t72 - t135 * t77) * pkin(4) + t3, -t94 * t72 + t106 + t68, t145 * t42 - t16 * t94 - t17 * t91 - t41 * t85, t72 * t91 + t77 * t85 + t1, t1 * t91 + t10 * t145 + t2 * t94 + t85 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, -t152, 0, -t103 * t117, t101 * t117, 0, 0, -t52, -t53, 0, -t27, t26, -t27, t145 * t81 - t52 * t94 - t53 * t91 - t80 * t85, -t26, t145 * t47 - t26 * t91 + t27 * t94 + t48 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t151, 0, 0, 0, 0, 0, -t53, t52, -t53, 0, -t52, t145 * t80 - t52 * t91 + t53 * t94 + t81 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -0.2e1 * t95, t89, 0, 0.2e1 * t85, 0.2e1 * t145 * t94 + 0.2e1 * t85 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, t72, t4, t3, t4 + 0.2e1 * t68, pkin(5) * t16 - qJ(6) * t17 - qJD(6) * t41, -t3 + 0.2e1 * t67 + 0.2e1 * t74, -pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t53, 0, -t27, t26, -t27, pkin(5) * t52 - qJ(6) * t53 - qJD(6) * t80, -t26, -pkin(5) * t27 - qJ(6) * t26 + qJD(6) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t52, -t53, 0, -t52, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t95, -t145, 0, t105 + t95, -pkin(5) * t145 + qJ(6) * t85 + qJD(6) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, qJ(6) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t16, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
