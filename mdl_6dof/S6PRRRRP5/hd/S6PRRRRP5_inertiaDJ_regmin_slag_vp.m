% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:26:03
% EndTime: 2019-03-09 00:26:11
% DurationCPUTime: 2.88s
% Computational Cost: add. (3243->341), mult. (9534->624), div. (0->0), fcn. (9405->12), ass. (0->169)
t103 = sin(qJ(4));
t201 = -0.4e1 * t103;
t107 = cos(qJ(4));
t108 = cos(qJ(3));
t169 = qJD(4) * t108;
t148 = t107 * t169;
t104 = sin(qJ(3));
t176 = qJD(3) * t104;
t98 = sin(pkin(7));
t200 = (t103 * t176 - t148) * t98;
t102 = sin(qJ(5));
t106 = cos(qJ(5));
t96 = t106 ^ 2;
t191 = t102 ^ 2 - t96;
t139 = t191 * qJD(5);
t199 = 0.2e1 * t98;
t97 = t107 ^ 2;
t198 = pkin(2) * t104;
t197 = pkin(10) * t102;
t196 = pkin(11) * t107;
t195 = -qJ(6) - pkin(11);
t100 = cos(pkin(7));
t189 = t104 * t98;
t61 = pkin(9) * t189 + (-pkin(2) * t108 - pkin(3)) * t100;
t70 = -t100 * t107 + t103 * t189;
t71 = t100 * t103 + t107 * t189;
t36 = t70 * pkin(4) - t71 * pkin(11) + t61;
t188 = t108 * t98;
t162 = pkin(9) * t188;
t62 = t162 + (pkin(10) + t198) * t100;
t131 = -pkin(3) * t108 - pkin(10) * t104;
t63 = (-pkin(2) + t131) * t98;
t194 = t103 * t63 + t107 * t62;
t38 = -pkin(11) * t188 + t194;
t13 = t102 * t36 + t106 * t38;
t129 = pkin(4) * t103 - t196;
t118 = t129 * qJD(4);
t167 = qJD(5) * t106;
t130 = -pkin(4) * t107 - pkin(11) * t103;
t82 = -pkin(3) + t130;
t193 = -t102 * t118 - t167 * t82;
t172 = qJD(4) * t103;
t145 = t102 * t172;
t192 = pkin(10) * t145 + t106 * t118;
t178 = t106 * t107;
t90 = pkin(10) * t178;
t60 = t102 * t82 + t90;
t95 = t103 ^ 2;
t190 = t95 - t97;
t109 = cos(qJ(2));
t99 = sin(pkin(6));
t187 = t109 * t99;
t112 = qJD(4) * t70;
t174 = qJD(3) * t108;
t150 = t107 * t174;
t134 = t98 * t150;
t111 = -t112 + t134;
t154 = t98 * t176;
t46 = t102 * t71 + t106 * t188;
t25 = qJD(5) * t46 - t102 * t154 - t106 * t111;
t186 = t25 * t102;
t185 = t25 * t106;
t184 = qJ(6) * t103;
t183 = qJD(2) * t99;
t182 = t103 * t106;
t105 = sin(qJ(2));
t181 = t104 * t105;
t180 = t104 * t109;
t179 = t105 * t108;
t177 = t108 * t109;
t175 = qJD(3) * t107;
t173 = qJD(4) * t102;
t171 = qJD(4) * t106;
t170 = qJD(4) * t107;
t168 = qJD(5) * t102;
t166 = qJD(5) * t107;
t165 = t106 * qJD(6);
t164 = -0.2e1 * pkin(3) * qJD(4);
t163 = -0.2e1 * pkin(4) * qJD(5);
t161 = t100 * t198;
t160 = t107 * t197;
t159 = pkin(10) * t170;
t158 = pkin(5) * t168;
t157 = t102 * t188;
t93 = t98 ^ 2;
t156 = t93 * t174;
t101 = cos(pkin(6));
t125 = t100 * t187 + t101 * t98;
t44 = t104 * t125 + t179 * t99;
t69 = t100 * t101 - t187 * t98;
t34 = t103 * t44 - t107 * t69;
t155 = t34 * t168;
t153 = t98 * t174;
t152 = t105 * t183;
t151 = t100 * t174;
t149 = t100 * t170;
t147 = t102 * t166;
t146 = t106 * t166;
t144 = t102 * t167;
t143 = t103 * t170;
t142 = t106 * t170;
t12 = -t102 * t38 + t106 * t36;
t141 = -t103 * t62 + t107 * t63;
t140 = qJD(5) * t195;
t138 = t190 * qJD(4);
t137 = 0.2e1 * t143;
t136 = t93 * t152;
t135 = t98 * t152;
t133 = t104 * t156;
t132 = t102 * t142;
t37 = pkin(4) * t188 - t141;
t47 = t106 * t71 - t157;
t8 = pkin(5) * t70 - qJ(6) * t47 + t12;
t9 = -qJ(6) * t46 + t13;
t128 = -t102 * t9 - t106 * t8;
t35 = t103 * t69 + t107 * t44;
t43 = -t108 * t125 + t181 * t99;
t17 = -t102 * t35 + t106 * t43;
t18 = t102 * t43 + t106 * t35;
t127 = -t102 * t18 - t106 * t17;
t126 = -t102 * t47 - t106 * t46;
t64 = (pkin(3) * t104 - pkin(10) * t108) * t98 * qJD(3);
t65 = -pkin(2) * t151 + pkin(9) * t154;
t20 = t103 * t65 + t107 * t64 - t170 * t62 - t172 * t63;
t30 = t101 * t153 + ((t100 * t177 - t181) * qJD(3) + (-t100 * t181 + t177) * qJD(2)) * t99;
t10 = qJD(4) * t35 + t30 * t103 - t107 * t135;
t124 = t10 * t102 + t167 * t34;
t123 = -t10 * t106 + t155;
t16 = -pkin(4) * t154 - t20;
t122 = t16 * t102 + t167 * t37;
t121 = -t16 * t106 + t168 * t37;
t45 = t71 * qJD(4) + t103 * t153;
t120 = t102 * t45 + t167 * t70;
t119 = -t106 * t45 + t168 * t70;
t19 = -t103 * t64 + t107 * t65 - t170 * t63 + t172 * t62;
t110 = t45 * pkin(4) + pkin(11) * t112 + (t161 + (pkin(9) - t196) * t188) * qJD(3);
t15 = pkin(11) * t154 - t19;
t3 = -t102 * t110 - t106 * t15 - t167 * t36 + t168 * t38;
t117 = -t104 * t172 + t150;
t115 = -t103 * t168 + t142;
t114 = t103 * t171 + t147;
t113 = t102 * t170 + t103 * t167;
t4 = -qJD(5) * t13 - t102 * t15 + t106 * t110;
t92 = -pkin(5) * t106 - pkin(4);
t84 = t195 * t106;
t83 = t195 * t102;
t78 = (pkin(5) * t102 + pkin(10)) * t103;
t77 = t106 * t82;
t68 = -t102 * qJD(6) + t106 * t140;
t67 = t102 * t140 + t165;
t66 = (t161 + t162) * qJD(3);
t59 = t77 - t160;
t55 = pkin(5) * t113 + t159;
t48 = -t102 * t184 + t60;
t42 = -qJ(6) * t182 + t77 + (-pkin(5) - t197) * t107;
t40 = -qJD(5) * t60 + t192;
t39 = pkin(10) * t114 + t193;
t29 = t101 * t154 + ((t100 * t180 + t179) * qJD(3) + (t100 * t179 + t180) * qJD(2)) * t99;
t27 = (-pkin(10) * qJD(4) - qJ(6) * qJD(5)) * t182 + (-qJD(6) * t103 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t107) * t102 - t193;
t26 = -qJD(5) * t157 + t102 * t111 - t106 * t154 + t167 * t71;
t24 = -t103 * t165 + (pkin(5) * t103 - qJ(6) * t178) * qJD(4) + (-t90 + (-t82 + t184) * t102) * qJD(5) + t192;
t21 = pkin(5) * t46 + t37;
t11 = -qJD(4) * t34 + t103 * t135 + t30 * t107;
t7 = pkin(5) * t26 + t16;
t6 = qJD(5) * t17 + t29 * t102 + t11 * t106;
t5 = -qJD(5) * t18 - t11 * t102 + t29 * t106;
t2 = -qJ(6) * t26 - qJD(6) * t46 - t3;
t1 = t45 * pkin(5) + t25 * qJ(6) - t47 * qJD(6) + t4;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t34 + 0.2e1 * t17 * t5 + 0.2e1 * t18 * t6; 0, 0, -t152, -t109 * t183, 0, 0, 0, 0, 0, -t100 * t29 - t108 * t136 + t154 * t69, -t100 * t30 + t104 * t136 + t153 * t69, 0, 0, 0, 0, 0, t29 * t70 + t43 * t45 + (t10 * t108 - t176 * t34) * t98, t43 * t149 + t29 * t71 + (t11 * t108 + t117 * t43 - t176 * t35) * t98, 0, 0, 0, 0, 0, t10 * t46 + t17 * t45 + t26 * t34 + t5 * t70, t10 * t47 - t18 * t45 - t25 * t34 - t6 * t70, t17 * t25 - t18 * t26 - t46 * t6 - t47 * t5, t1 * t17 + t10 * t21 + t18 * t2 + t34 * t7 + t5 * t8 + t6 * t9; 0, 0, 0, 0, 0.2e1 * t133, 0.2e1 * (-t104 ^ 2 + t108 ^ 2) * t93 * qJD(3), t151 * t199, -0.2e1 * t100 * t154, 0, -0.2e1 * pkin(2) * t176 * t93 - 0.2e1 * t100 * t66, -0.2e1 * pkin(2) * t156 + 0.2e1 * t100 * t65, 0.2e1 * t71 * t111, -0.2e1 * t111 * t70 - 0.2e1 * t71 * t45 (-t100 * t148 - t117 * t188 + t176 * t71) * t199 (t108 * t45 - t176 * t70) * t199, -0.2e1 * t133, 0.2e1 * t61 * t45 + 0.2e1 * t66 * t70 + 0.2e1 * (-t20 * t108 + t141 * t176) * t98, 0.2e1 * t61 * t149 + 0.2e1 * t66 * t71 + 0.2e1 * (-t19 * t108 + t117 * t61 - t176 * t194) * t98, -0.2e1 * t47 * t25, 0.2e1 * t25 * t46 - 0.2e1 * t26 * t47, -0.2e1 * t25 * t70 + 0.2e1 * t45 * t47, -0.2e1 * t26 * t70 - 0.2e1 * t45 * t46, 0.2e1 * t70 * t45, 0.2e1 * t12 * t45 + 0.2e1 * t16 * t46 + 0.2e1 * t26 * t37 + 0.2e1 * t4 * t70, -0.2e1 * t13 * t45 + 0.2e1 * t16 * t47 - 0.2e1 * t25 * t37 + 0.2e1 * t3 * t70, -0.2e1 * t1 * t47 - 0.2e1 * t2 * t46 + 0.2e1 * t25 * t8 - 0.2e1 * t26 * t9, 0.2e1 * t1 * t8 + 0.2e1 * t2 * t9 + 0.2e1 * t21 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t30, 0, 0, 0, 0, 0, -t107 * t29 + t172 * t43, t103 * t29 + t170 * t43, 0, 0, 0, 0, 0 (t173 * t34 - t5) * t107 + (qJD(4) * t17 + t124) * t103 (t171 * t34 + t6) * t107 + (-qJD(4) * t18 - t123) * t103, t127 * t170 + (-t102 * t6 - t106 * t5 + (t102 * t17 - t106 * t18) * qJD(5)) * t103, t10 * t78 + t17 * t24 + t18 * t27 + t34 * t55 + t42 * t5 + t48 * t6; 0, 0, 0, 0, 0, 0, t153, -t154, 0, -t66, t65, t103 * t134 + (-t103 * t70 + t71 * t107) * qJD(4), t97 * t153 - t103 * t45 + (-t71 * t103 - 0.2e1 * t107 * t70) * qJD(4), t200 (t103 * t169 + t104 * t175) * t98, 0, -pkin(3) * t45 - pkin(10) * t200 - t66 * t107 + t172 * t61, t66 * t103 + t131 * t98 * t175 + (-pkin(10) * t103 * t188 + pkin(3) * t70 + t61 * t107) * qJD(4), t47 * t142 + (-t168 * t47 - t185) * t103, t126 * t170 + (t186 - t106 * t26 + (t102 * t46 - t106 * t47) * qJD(5)) * t103 (t171 * t70 + t25) * t107 + (qJD(4) * t47 - t119) * t103 (-t173 * t70 + t26) * t107 + (-qJD(4) * t46 - t120) * t103, -t107 * t45 + t172 * t70, t40 * t70 + t59 * t45 + (-t4 + (pkin(10) * t46 + t102 * t37) * qJD(4)) * t107 + (pkin(10) * t26 + qJD(4) * t12 + t122) * t103, t39 * t70 - t60 * t45 + (-t3 + (pkin(10) * t47 + t106 * t37) * qJD(4)) * t107 + (-pkin(10) * t25 - qJD(4) * t13 - t121) * t103, -t24 * t47 + t42 * t25 - t48 * t26 - t27 * t46 + t128 * t170 + (-t1 * t106 - t102 * t2 + (t102 * t8 - t106 * t9) * qJD(5)) * t103, t1 * t42 + t2 * t48 + t21 * t55 + t24 * t8 + t27 * t9 + t7 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, -0.2e1 * t138, 0, 0, 0, t103 * t164, t107 * t164, 0.2e1 * t143 * t96 - 0.2e1 * t144 * t95, t132 * t201 + 0.2e1 * t139 * t95, 0.2e1 * t103 * t147 + 0.2e1 * t171 * t190, -0.2e1 * t102 * t138 + 0.2e1 * t103 * t146, -0.2e1 * t143, 0.2e1 * t59 * t172 - 0.2e1 * t40 * t107 + 0.2e1 * (t102 * t137 + t167 * t95) * pkin(10), -0.2e1 * t60 * t172 - 0.2e1 * t39 * t107 + 0.2e1 * (t106 * t137 - t168 * t95) * pkin(10), 0.2e1 * (-t102 * t48 - t106 * t42) * t170 + 0.2e1 * (-t102 * t27 - t106 * t24 + (t102 * t42 - t106 * t48) * qJD(5)) * t103, 0.2e1 * t24 * t42 + 0.2e1 * t27 * t48 + 0.2e1 * t55 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, t123, t124, qJD(5) * t127 - t5 * t102 + t6 * t106, pkin(5) * t155 + t10 * t92 + t17 * t68 + t18 * t67 + t5 * t83 - t6 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t45, t154, t20, t19, t167 * t47 - t186, qJD(5) * t126 - t102 * t26 - t185, t120, -t119, 0, -pkin(4) * t26 - pkin(11) * t120 + t121, pkin(4) * t25 + pkin(11) * t119 + t122, qJD(5) * t128 - t1 * t102 + t2 * t106 + t83 * t25 + t84 * t26 - t67 * t46 - t68 * t47, t1 * t83 + t158 * t21 - t2 * t84 + t67 * t9 + t68 * t8 + t7 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, -t172, 0, -t159, pkin(10) * t172, -t103 * t139 + t132, t144 * t201 - t170 * t191, t145 - t146, t114, 0 (pkin(11) * t178 + (-pkin(4) * t106 + t197) * t103) * qJD(5) + (t102 * t130 - t90) * qJD(4) (pkin(10) * t182 + t102 * t129) * qJD(5) + (t106 * t130 + t160) * qJD(4) (-t83 * t170 - t103 * t68 + t27 + (t103 * t84 - t42) * qJD(5)) * t106 + (t84 * t170 - t103 * t67 - t24 + (t103 * t83 - t48) * qJD(5)) * t102, t158 * t78 + t24 * t83 - t27 * t84 + t42 * t68 + t48 * t67 + t55 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t144, -0.2e1 * t139, 0, 0, 0, t102 * t163, t106 * t163, -0.2e1 * t68 * t102 + 0.2e1 * t67 * t106 + 0.2e1 * (t102 * t84 - t106 * t83) * qJD(5), 0.2e1 * t158 * t92 - 0.2e1 * t67 * t84 + 0.2e1 * t68 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, t45, t4, t3, pkin(5) * t25, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t113, t172, t40, t39, -t115 * pkin(5), t24 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, -t168, 0, -pkin(11) * t167, pkin(11) * t168, -pkin(5) * t167, t68 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t14;