% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:52
% EndTime: 2019-12-31 22:29:59
% DurationCPUTime: 2.01s
% Computational Cost: add. (2345->214), mult. (5965->401), div. (0->0), fcn. (5332->8), ass. (0->136)
t103 = cos(qJ(2));
t143 = t103 * qJD(2);
t98 = sin(qJ(3));
t131 = t98 * t143;
t102 = cos(qJ(3));
t147 = qJD(3) * t102;
t99 = sin(qJ(2));
t172 = t99 * t147 + t131;
t171 = -0.4e1 * t99;
t101 = cos(qJ(4));
t97 = sin(qJ(4));
t67 = t101 * t98 + t102 * t97;
t54 = t67 * t99;
t153 = t102 * t99;
t167 = pkin(6) * t98;
t112 = -pkin(2) * t103 - t99 * pkin(7);
t75 = -pkin(1) + t112;
t65 = t102 * t75;
t45 = -pkin(8) * t153 + t65 + (-pkin(3) - t167) * t103;
t149 = t102 * t103;
t85 = pkin(6) * t149;
t158 = t98 * t75 + t85;
t163 = t98 * t99;
t50 = -pkin(8) * t163 + t158;
t47 = t101 * t50;
t162 = t97 * t45 + t47;
t168 = pkin(7) + pkin(8);
t78 = t168 * t98;
t79 = t168 * t102;
t159 = t101 * t79 - t97 * t78;
t124 = t102 * t143;
t152 = qJD(3) * t98;
t170 = -t99 * t152 + t124;
t94 = t102 ^ 2;
t157 = t98 ^ 2 - t94;
t117 = t157 * qJD(3);
t169 = qJD(3) + qJD(4);
t146 = qJD(3) * t103;
t129 = t98 * t146;
t148 = qJD(2) * t102;
t106 = t99 * t148 + t129;
t113 = pkin(2) * t99 - pkin(7) * t103;
t73 = t113 * qJD(2);
t30 = pkin(6) * t106 - t75 * t147 - t98 * t73;
t166 = pkin(3) * t101;
t165 = pkin(4) * t103;
t96 = sin(qJ(5));
t164 = t96 * t97;
t90 = t99 * qJD(2);
t133 = t98 * t90;
t160 = pkin(6) * t133 + t102 * t73;
t74 = pkin(3) * t163 + t99 * pkin(6);
t93 = t99 ^ 2;
t156 = -t103 ^ 2 + t93;
t100 = cos(qJ(5));
t18 = -pkin(9) * t54 + t162;
t155 = t100 * t18;
t154 = t100 * t97;
t151 = qJD(4) * t97;
t150 = qJD(5) * t96;
t145 = qJD(4) * t101;
t144 = qJD(5) * t100;
t142 = -0.2e1 * pkin(1) * qJD(2);
t141 = -0.2e1 * pkin(2) * qJD(3);
t89 = pkin(6) * t143;
t51 = t172 * pkin(3) + t89;
t140 = pkin(3) * t152;
t139 = pkin(4) * t90;
t138 = pkin(3) * t151;
t137 = pkin(4) * t150;
t135 = pkin(3) * t145;
t134 = pkin(4) * t144;
t108 = t101 * t102 - t97 * t98;
t26 = t108 * t143 - t169 * t54;
t22 = (pkin(3) * t99 - pkin(8) * t149) * qJD(2) + (-t85 + (pkin(8) * t99 - t75) * t98) * qJD(3) + t160;
t25 = -pkin(8) * t172 - t30;
t121 = t101 * t22 - t97 * t25;
t9 = -t162 * qJD(4) + t121;
t6 = -pkin(9) * t26 + t139 + t9;
t27 = -t151 * t163 + (t169 * t153 + t131) * t101 + t170 * t97;
t8 = -t101 * t25 - t45 * t145 + t50 * t151 - t97 * t22;
t7 = -pkin(9) * t27 - t8;
t132 = t100 * t6 - t96 * t7;
t128 = t98 * t147;
t127 = t99 * t143;
t88 = -pkin(3) * t102 - pkin(2);
t14 = t18 * t150;
t126 = -t96 * t6 + t14;
t125 = qJD(3) * t168;
t123 = t102 * t146;
t120 = t101 * t45 - t97 * t50;
t55 = t108 * t99;
t17 = -t55 * pkin(9) + t120 - t165;
t122 = -t17 + t165;
t119 = -t101 * t78 - t79 * t97;
t87 = pkin(4) + t166;
t118 = qJD(5) * (-pkin(4) - t87);
t116 = t156 * qJD(2);
t115 = 0.2e1 * t127;
t114 = t98 * t124;
t111 = t17 * t96 + t155;
t34 = -pkin(9) * t67 + t119;
t35 = pkin(9) * t108 + t159;
t110 = t100 * t35 + t34 * t96;
t109 = t100 * t34 - t35 * t96;
t33 = t100 * t55 - t54 * t96;
t32 = t100 * t54 + t55 * t96;
t44 = t100 * t67 + t108 * t96;
t43 = -t100 * t108 + t67 * t96;
t71 = t98 * t125;
t72 = t102 * t125;
t28 = t101 * t71 + t78 * t145 + t79 * t151 + t97 * t72;
t2 = -qJD(5) * t111 + t132;
t29 = -qJD(4) * t159 - t101 * t72 + t97 * t71;
t105 = (t97 * t150 + (-t100 * t101 + t164) * qJD(4)) * pkin(3);
t104 = (-t97 * t144 + (-t101 * t96 - t154) * qJD(4)) * pkin(3);
t83 = -0.2e1 * t127;
t53 = -pkin(4) * t108 + t88;
t49 = t169 * t67;
t48 = t169 * t108;
t46 = pkin(4) * t54 + t74;
t39 = -t150 * t87 + t104;
t38 = -t144 * t87 + t105;
t36 = pkin(4) * t49 + t140;
t31 = -qJD(3) * t158 + t160;
t19 = pkin(4) * t27 + t51;
t16 = -pkin(9) * t48 + t29;
t15 = -pkin(9) * t49 - t28;
t13 = qJD(5) * t44 + t100 * t49 + t96 * t48;
t12 = -qJD(5) * t43 + t100 * t48 - t96 * t49;
t11 = qJD(5) * t33 + t100 * t27 + t96 * t26;
t10 = -qJD(5) * t32 + t100 * t26 - t96 * t27;
t5 = -qJD(5) * t110 + t100 * t16 - t96 * t15;
t4 = -qJD(5) * t109 - t100 * t15 - t96 * t16;
t1 = (-qJD(5) * t17 - t7) * t100 + t126;
t3 = [0, 0, 0, t115, -0.2e1 * t116, 0, 0, 0, t99 * t142, t103 * t142, 0.2e1 * t127 * t94 - 0.2e1 * t128 * t93, t114 * t171 + 0.2e1 * t117 * t93, 0.2e1 * t129 * t99 + 0.2e1 * t148 * t156, -0.2e1 * t116 * t98 + 0.2e1 * t123 * t99, t83, 0.2e1 * t65 * t90 - 0.2e1 * t31 * t103 + 0.2e1 * (t127 * t98 + t147 * t93) * pkin(6), -0.2e1 * t30 * t103 - 0.2e1 * t158 * t90 + 0.2e1 * (t102 * t115 - t152 * t93) * pkin(6), 0.2e1 * t55 * t26, -0.2e1 * t26 * t54 - 0.2e1 * t27 * t55, -0.2e1 * t103 * t26 + 0.2e1 * t55 * t90, 0.2e1 * t103 * t27 - 0.2e1 * t54 * t90, t83, -0.2e1 * t103 * t9 + 0.2e1 * t120 * t90 + 0.2e1 * t27 * t74 + 0.2e1 * t51 * t54, -0.2e1 * t103 * t8 - 0.2e1 * t162 * t90 + 0.2e1 * t26 * t74 + 0.2e1 * t51 * t55, 0.2e1 * t33 * t10, -0.2e1 * t10 * t32 - 0.2e1 * t11 * t33, -0.2e1 * t10 * t103 + 0.2e1 * t33 * t90, 0.2e1 * t103 * t11 - 0.2e1 * t32 * t90, t83, -0.2e1 * t2 * t103 + 0.2e1 * (t100 * t17 - t18 * t96) * t90 + 0.2e1 * t19 * t32 + 0.2e1 * t46 * t11, -0.2e1 * t1 * t103 + 0.2e1 * t10 * t46 - 0.2e1 * t111 * t90 + 0.2e1 * t19 * t33; 0, 0, 0, 0, 0, t143, -t90, 0, -t89, pkin(6) * t90, -t117 * t99 + t114, t128 * t171 - t143 * t157, -t123 + t133, t106, 0, (pkin(6) * t163 - t102 * t113) * qJD(3) + (-pkin(7) * t163 + (-pkin(2) * t98 - pkin(6) * t102) * t103) * qJD(2), (pkin(6) * t153 + t113 * t98) * qJD(3) + (t102 * t112 + t103 * t167) * qJD(2), t26 * t67 + t48 * t55, t108 * t26 - t27 * t67 - t48 * t54 - t49 * t55, -t103 * t48 + t67 * t90, t103 * t49 + t108 * t90, 0, -t103 * t29 - t108 * t51 + t119 * t90 + t140 * t54 + t27 * t88 + t49 * t74, -t103 * t28 + t140 * t55 - t159 * t90 + t26 * t88 + t48 * t74 + t51 * t67, t10 * t44 + t12 * t33, -t10 * t43 - t11 * t44 - t12 * t32 - t13 * t33, -t103 * t12 + t44 * t90, t103 * t13 - t43 * t90, 0, -t103 * t5 + t109 * t90 + t11 * t53 + t13 * t46 + t19 * t43 + t32 * t36, t10 * t53 - t103 * t4 - t110 * t90 + t12 * t46 + t19 * t44 + t33 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t128, -0.2e1 * t117, 0, 0, 0, t98 * t141, t102 * t141, 0.2e1 * t67 * t48, 0.2e1 * t108 * t48 - 0.2e1 * t49 * t67, 0, 0, 0, -0.2e1 * t108 * t140 + 0.2e1 * t49 * t88, 0.2e1 * t140 * t67 + 0.2e1 * t48 * t88, 0.2e1 * t44 * t12, -0.2e1 * t12 * t43 - 0.2e1 * t13 * t44, 0, 0, 0, 0.2e1 * t13 * t53 + 0.2e1 * t36 * t43, 0.2e1 * t12 * t53 + 0.2e1 * t36 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, -t172, t90, t31, t30, 0, 0, t26, -t27, t90, t90 * t166 + (-t47 + (pkin(3) * t103 - t45) * t97) * qJD(4) + t121, (t103 * t145 - t90 * t97) * pkin(3) + t8, 0, 0, t10, -t11, t90, -t39 * t103 + (-pkin(3) * t164 + t100 * t87) * t90 + t2, -t38 * t103 - (pkin(3) * t154 + t87 * t96) * t90 - t100 * t7 - t17 * t144 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, -t152, 0, -pkin(7) * t147, pkin(7) * t152, 0, 0, t48, -t49, 0, t29, t28, 0, 0, t12, -t13, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t138, -0.2e1 * t135, 0, 0, 0, 0, 0, 0.2e1 * t39, 0.2e1 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t27, t90, t9, t8, 0, 0, t10, -t11, t90, t100 * t139 + (t122 * t96 - t155) * qJD(5) + t132, t14 + (-t6 - t139) * t96 + (qJD(5) * t122 - t7) * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t49, 0, t29, t28, 0, 0, t12, -t13, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, -t135, 0, 0, 0, 0, 0, t118 * t96 + t104, t100 * t118 + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t137, -0.2e1 * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t90, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t13, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, -t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
