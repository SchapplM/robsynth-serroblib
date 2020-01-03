% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:22
% EndTime: 2019-12-31 19:36:29
% DurationCPUTime: 2.17s
% Computational Cost: add. (2961->281), mult. (7530->375), div. (0->0), fcn. (5210->6), ass. (0->160)
t104 = sin(pkin(8));
t105 = cos(pkin(8));
t107 = sin(qJ(2));
t109 = cos(qJ(2));
t84 = t104 * t109 + t105 * t107;
t116 = qJD(1) * t84;
t192 = qJD(5) + t116;
t106 = sin(qJ(5));
t108 = cos(qJ(5));
t160 = t105 * t109;
t143 = qJD(1) * t160;
t151 = qJD(1) * t107;
t73 = t104 * t151 - t143;
t50 = qJD(2) * t106 - t108 * t73;
t138 = t192 * t50;
t148 = qJD(5) * t108;
t149 = qJD(5) * t106;
t75 = t84 * qJD(2);
t64 = qJD(1) * t75;
t31 = qJD(2) * t149 - t106 * t64 - t73 * t148;
t193 = t31 - t138;
t147 = qJD(1) * qJD(2);
t141 = t109 * t147;
t142 = t107 * t147;
t91 = t104 * t142;
t65 = t105 * t141 - t91;
t95 = pkin(2) * t142;
t140 = -qJ(4) * t65 + t95;
t119 = -qJD(4) * t116 + t140;
t180 = pkin(3) + pkin(7);
t11 = t180 * t64 + t119;
t166 = -qJ(3) - pkin(6);
t139 = qJD(2) * t166;
t68 = qJD(3) * t109 + t107 * t139;
t59 = t68 * qJD(1);
t69 = -qJD(3) * t107 + t109 * t139;
t60 = t69 * qJD(1);
t27 = t104 * t59 - t105 * t60;
t15 = pkin(4) * t65 + t27;
t99 = -pkin(2) * t109 - pkin(1);
t162 = qJD(1) * t99;
t88 = qJD(3) + t162;
t114 = -qJ(4) * t116 + t88;
t17 = t180 * t73 + t114;
t90 = t166 * t109;
t87 = qJD(1) * t90;
t79 = t104 * t87;
t89 = t166 * t107;
t86 = qJD(1) * t89;
t82 = qJD(2) * pkin(2) + t86;
t43 = t105 * t82 + t79;
t130 = qJD(4) - t43;
t177 = pkin(4) * t116;
t18 = -t180 * qJD(2) + t130 + t177;
t6 = t106 * t18 + t108 * t17;
t2 = -qJD(5) * t6 - t106 * t11 + t108 * t15;
t184 = t192 * t6 + t2;
t122 = t106 * t17 - t108 * t18;
t1 = -t122 * qJD(5) + t106 * t15 + t108 * t11;
t129 = t122 * t192 + t1;
t185 = t108 * t192;
t52 = qJD(2) * t108 + t106 * t73;
t194 = t52 * t185;
t136 = t106 * t192;
t120 = t108 * t65 - t136 * t192;
t191 = -0.2e1 * t147;
t181 = t73 ^ 2;
t72 = t116 ^ 2;
t189 = -t181 - t72;
t188 = -t181 + t72;
t45 = t65 * t84;
t150 = qJD(2) * t107;
t78 = qJD(2) * t160 - t104 * t150;
t187 = t116 * t78 + t45;
t186 = 0.2e1 * t116;
t47 = t105 * t86 + t79;
t155 = -qJD(4) + t47;
t178 = pkin(4) * t73;
t165 = t105 * t87;
t44 = t104 * t82 - t165;
t40 = -qJD(2) * qJ(4) - t44;
t23 = -t40 - t178;
t98 = -pkin(2) * t105 - pkin(3);
t94 = -pkin(7) + t98;
t183 = t192 * t23 + t65 * t94;
t182 = (-t73 + t143) * qJD(2) - t91;
t179 = pkin(3) * t64;
t176 = pkin(2) * t107;
t48 = -t104 * t90 - t105 * t89;
t175 = t27 * t48;
t33 = pkin(3) * t73 + t114;
t174 = t33 * t116;
t173 = t50 * t73;
t172 = t52 * t50;
t171 = t52 * t73;
t83 = t104 * t107 - t160;
t170 = t65 * t83;
t168 = t73 * t116;
t28 = t104 * t60 + t105 * t59;
t57 = t108 * t64;
t32 = t52 * qJD(5) - t57;
t164 = t106 * t32;
t163 = t108 * t31;
t161 = qJD(2) * t75;
t111 = qJD(1) ^ 2;
t159 = t109 * t111;
t110 = qJD(2) ^ 2;
t158 = t110 * t107;
t157 = t110 * t109;
t156 = t78 * qJD(2);
t154 = t177 - t155;
t152 = t107 ^ 2 - t109 ^ 2;
t101 = pkin(2) * t150;
t100 = pkin(2) * t151;
t146 = t83 * t149;
t145 = t83 * t148;
t144 = t107 * t159;
t35 = t104 * t68 - t105 * t69;
t46 = t104 * t86 - t165;
t137 = qJ(4) * t73 + t100;
t133 = pkin(1) * t191;
t132 = t107 * t141;
t24 = -qJD(2) * qJD(4) - t28;
t13 = -pkin(4) * t64 - t24;
t131 = -qJD(5) * t192 * t94 + t13;
t127 = t106 * t122 + t108 * t6;
t126 = t13 * t83 + t23 * t75;
t125 = t64 * t83 + t73 * t75;
t124 = t192 * t75 + t170;
t36 = t104 * t69 + t105 * t68;
t49 = t104 * t89 - t105 * t90;
t121 = -qJ(4) * t84 + t99;
t26 = t180 * t83 + t121;
t37 = pkin(4) * t84 + t48;
t10 = t106 * t37 + t108 * t26;
t9 = -t106 * t26 + t108 * t37;
t118 = -qJD(2) * t46 + t27;
t117 = -qJ(4) * t78 - qJD(4) * t84 + t101;
t115 = -t106 * t65 - t185 * t192;
t41 = -t91 + (t73 + t143) * qJD(2);
t113 = -t116 * t75 - t64 * t84 - t73 * t78 - t170;
t112 = t116 * t35 + t27 * t84 - t36 * t73 + t48 * t65 - t49 * t64;
t96 = pkin(2) * t104 + qJ(4);
t42 = pkin(3) * t83 + t121;
t39 = -qJD(2) * pkin(3) + t130;
t38 = -pkin(4) * t83 + t49;
t34 = pkin(3) * t116 + t137;
t29 = t46 - t178;
t25 = t108 * t32;
t22 = pkin(3) * t75 + t117;
t21 = -pkin(4) * t75 + t36;
t20 = pkin(4) * t78 + t35;
t19 = t116 * t180 + t137;
t16 = t119 + t179;
t12 = t180 * t75 + t117;
t8 = t106 * t29 + t108 * t19;
t7 = -t106 * t19 + t108 * t29;
t4 = -t10 * qJD(5) - t106 * t12 + t108 * t20;
t3 = t9 * qJD(5) + t106 * t20 + t108 * t12;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t132, t152 * t191, t157, -0.2e1 * t132, -t158, 0, -pkin(6) * t157 + t107 * t133, pkin(6) * t158 + t109 * t133, 0, 0, t187, t113, t156, t125, -t161, 0, t64 * t99 + t75 * t88 + (-t35 + (qJD(1) * t83 + t73) * t176) * qJD(2), t65 * t99 + t78 * t88 + (t186 * t176 - t36) * qJD(2), -t28 * t83 - t43 * t78 - t44 * t75 + t112, t175 + t28 * t49 - t35 * t43 + t36 * t44 + (t88 + t162) * t101, 0, -t156, t161, t187, t113, t125, t24 * t83 + t39 * t78 + t40 * t75 + t112, qJD(2) * t35 - t16 * t83 - t22 * t73 - t33 * t75 - t42 * t64, qJD(2) * t36 - t116 * t22 - t16 * t84 - t33 * t78 - t42 * t65, t16 * t42 + t22 * t33 - t24 * t49 + t35 * t39 - t36 * t40 + t175, t52 * t145 + (-t31 * t83 + t52 * t75) * t106, (-t106 * t50 + t108 * t52) * t75 + (-t164 - t163 + (-t106 * t52 - t108 * t50) * qJD(5)) * t83, t106 * t124 + t145 * t192 - t31 * t84 + t52 * t78, t50 * t146 + (-t32 * t83 - t50 * t75) * t108, t108 * t124 - t146 * t192 - t32 * t84 - t50 * t78, t192 * t78 + t45, -t108 * t126 - t122 * t78 + t146 * t23 + t192 * t4 + t2 * t84 + t21 * t50 + t32 * t38 + t65 * t9, -t1 * t84 - t10 * t65 + t106 * t126 + t145 * t23 - t192 * t3 + t21 * t52 - t31 * t38 - t6 * t78, -t10 * t32 - t3 * t50 + t31 * t9 - t4 * t52 + t127 * t75 + (t1 * t108 - t106 * t2 + (-t106 * t6 + t108 * t122) * qJD(5)) * t83, t1 * t10 - t122 * t4 + t13 * t38 + t2 * t9 + t21 * t23 + t3 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t152 * t111, 0, t144, 0, 0, t111 * pkin(1) * t107, pkin(1) * t159, 0, 0, t168, t188, t41, -t168, 0, 0, -t73 * t100 - t116 * t88 - t118, qJD(2) * t47 - t100 * t116 + t73 * t88 - t28, (t44 - t46) * t116 + (-t43 + t47) * t73 + (-t104 * t64 - t105 * t65) * pkin(2), t43 * t46 - t44 * t47 + (t104 * t28 - t105 * t27 - t88 * t151) * pkin(2), 0, -t41, 0, t168, t188, -t168, -t64 * t96 + t65 * t98 + (-t40 - t46) * t116 + (t39 + t155) * t73, t34 * t73 + t118 + t174, -t33 * t73 + t34 * t116 + (0.2e1 * qJD(4) - t47) * qJD(2) + t28, t155 * t40 - t24 * t96 + t27 * t98 - t33 * t34 - t39 * t46, -t136 * t52 - t163, -t25 - t194 + (t31 + t138) * t106, t120 + t171, t185 * t50 + t164, t115 - t173, t192 * t73, t131 * t106 + t183 * t108 - t122 * t73 + t154 * t50 - t192 * t7 + t32 * t96, -t183 * t106 + t131 * t108 + t154 * t52 + t192 * t8 - t31 * t96 - t6 * t73, t50 * t8 + t52 * t7 + (t31 * t94 - t6 * t116 - t2 + (-t50 * t94 - t6) * qJD(5)) * t108 + (-t32 * t94 - t122 * t116 - t1 + (t52 * t94 - t122) * qJD(5)) * t106, t13 * t96 + t122 * t7 - t6 * t8 + t154 * t23 + (qJD(5) * t127 + t1 * t106 + t108 * t2) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186 * qJD(2), t182, t189, t116 * t43 + t44 * t73 + t95, 0, 0, 0, 0, 0, 0, t189, -0.2e1 * qJD(2) * t116, -t182, t179 - t40 * t73 + (-qJD(4) - t39) * t116 + t140, 0, 0, 0, 0, 0, 0, t115 + t173, t171 - t120, -t193 * t106 + t194 - t25, -t184 * t106 + t129 * t108 + t23 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t168, -t72 - t110, qJD(2) * t40 + t174 + t27, 0, 0, 0, 0, 0, 0, -qJD(2) * t50 + t120, -qJD(2) * t52 + t115, t193 * t108 + (t192 * t52 - t32) * t106, -qJD(2) * t23 + t129 * t106 + t184 * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, -t50 ^ 2 + t52 ^ 2, -t193, -t172, t57 + (-qJD(5) + t192) * t52, t65, -t23 * t52 + t184, t23 * t50 - t129, 0, 0;];
tauc_reg = t5;
