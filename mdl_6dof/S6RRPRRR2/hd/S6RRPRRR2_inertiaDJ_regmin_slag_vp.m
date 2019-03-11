% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:09
% EndTime: 2019-03-09 13:19:16
% DurationCPUTime: 2.10s
% Computational Cost: add. (4653->230), mult. (10318->385), div. (0->0), fcn. (10880->10), ass. (0->172)
t129 = sin(qJ(5));
t132 = cos(qJ(6));
t128 = sin(qJ(6));
t133 = cos(qJ(5));
t186 = t128 * t133;
t107 = t132 * t129 + t186;
t127 = sin(pkin(11));
t131 = sin(qJ(2));
t134 = cos(qJ(2));
t188 = cos(pkin(11));
t101 = t127 * t134 + t188 * t131;
t130 = sin(qJ(4));
t158 = t188 * t134;
t142 = t127 * t131 - t158;
t211 = cos(qJ(4));
t74 = t211 * t101 - t130 * t142;
t46 = t107 * t74;
t155 = t188 * pkin(2) + pkin(3);
t210 = pkin(2) * t127;
t217 = t130 * t210 - t211 * t155;
t123 = qJD(5) * t133;
t169 = t74 * t123;
t179 = t131 * qJD(2);
t100 = qJD(2) * t158 - t127 * t179;
t137 = t211 * t142;
t140 = qJD(2) * t101;
t183 = qJD(4) * t130;
t48 = -qJD(4) * t137 + t211 * t100 - t101 * t183 - t130 * t140;
t195 = t129 * t48;
t216 = t169 + t195;
t126 = t133 ^ 2;
t184 = t129 ^ 2 - t126;
t156 = t184 * qJD(5);
t177 = qJD(5) + qJD(6);
t215 = -pkin(9) - pkin(10);
t49 = t74 * qJD(4) + t130 * t100 + t211 * t140;
t214 = t49 * pkin(5);
t73 = t130 * t101 + t137;
t213 = t73 * pkin(5);
t136 = t130 * t155 + t211 * t210;
t94 = pkin(9) + t136;
t212 = -pkin(10) - t94;
t182 = qJD(5) * t129;
t206 = -qJ(3) - pkin(7);
t159 = qJD(2) * t206;
t97 = t134 * qJD(3) + t131 * t159;
t98 = -t131 * qJD(3) + t134 * t159;
t63 = t127 * t98 + t188 * t97;
t135 = -pkin(8) * t140 + t63;
t62 = -t127 * t97 + t188 * t98;
t139 = -t100 * pkin(8) + t62;
t163 = qJD(4) * t211;
t111 = t206 * t131;
t112 = t206 * t134;
t77 = t188 * t111 + t127 * t112;
t64 = -t101 * pkin(8) + t77;
t78 = t127 * t111 - t188 * t112;
t65 = -t142 * pkin(8) + t78;
t19 = -t130 * t139 - t211 * t135 - t64 * t163 + t65 * t183;
t122 = pkin(2) * t179;
t82 = pkin(3) * t140 + t122;
t24 = t49 * pkin(4) - t48 * pkin(9) + t82;
t168 = -t134 * pkin(2) - pkin(1);
t87 = t142 * pkin(3) + t168;
t42 = t73 * pkin(4) - t74 * pkin(9) + t87;
t44 = t130 * t64 + t211 * t65;
t7 = -t42 * t123 - t129 * t24 + t133 * t19 + t44 * t182;
t6 = -pkin(10) * t216 - t7;
t209 = t132 * t6;
t208 = t133 * pkin(5);
t207 = t74 * t48;
t187 = t128 * t129;
t106 = -t132 * t133 + t187;
t20 = t130 * t135 - t211 * t139 + t65 * t163 + t64 * t183;
t12 = pkin(5) * t216 + t20;
t194 = t129 * t74;
t43 = t130 * t65 - t211 * t64;
t32 = pkin(5) * t194 + t43;
t76 = t177 * t107;
t205 = t12 * t106 + t32 * t76;
t180 = qJD(6) * t132;
t75 = -t132 * t123 - t133 * t180 + t177 * t187;
t204 = t12 * t107 - t32 * t75;
t38 = t43 * t123;
t203 = t20 * t129 + t38;
t41 = t133 * t44;
t202 = t129 * t42 + t41;
t173 = pkin(5) * t182;
t91 = t136 * qJD(4);
t81 = t91 + t173;
t93 = -pkin(4) + t217;
t86 = t93 - t208;
t201 = t81 * t106 + t86 * t76;
t200 = t81 * t107 - t86 * t75;
t199 = t93 * t123 + t91 * t129;
t121 = -pkin(4) - t208;
t198 = t106 * t173 + t121 * t76;
t197 = t107 * t173 - t121 * t75;
t16 = -pkin(10) * t194 + t202;
t196 = t128 * t16;
t90 = t217 * qJD(4);
t193 = t129 * t90;
t192 = t132 * t16;
t191 = t133 * t48;
t190 = t133 * t74;
t189 = t133 * t90;
t185 = t129 * t133;
t181 = qJD(6) * t128;
t178 = t134 * qJD(2);
t176 = -0.2e1 * pkin(1) * qJD(2);
t175 = pkin(4) * t182;
t174 = pkin(4) * t123;
t172 = pkin(5) * t181;
t171 = pkin(5) * t180;
t170 = t74 * t182;
t37 = t43 * t182;
t40 = t133 * t42;
t15 = -pkin(10) * t190 - t129 * t44 + t213 + t40;
t167 = -t15 - t213;
t161 = t129 * t19 + t133 * t24;
t5 = -pkin(10) * t191 + t214 + (-t41 + (pkin(10) * t74 - t42) * t129) * qJD(5) + t161;
t166 = -t128 * t6 + t132 * t5;
t165 = qJD(5) * t215;
t164 = t129 * t123;
t162 = qJD(5) * t212;
t160 = -t91 * t133 + t93 * t182;
t157 = -0.4e1 * t74 * t185;
t154 = t20 * t74 + t43 * t48;
t153 = t48 * t73 + t74 * t49;
t152 = t73 * t94 - t74 * t93;
t30 = t107 * t49 - t75 * t73;
t151 = -t132 * t15 + t196;
t150 = t128 * t15 + t192;
t79 = t212 * t129;
t124 = t133 * pkin(10);
t80 = t133 * t94 + t124;
t149 = t128 * t80 - t132 * t79;
t148 = t128 * t79 + t132 * t80;
t113 = t215 * t129;
t114 = t133 * pkin(9) + t124;
t147 = t132 * t113 - t128 * t114;
t146 = t128 * t113 + t132 * t114;
t144 = t170 - t191;
t34 = t73 * t123 + t129 * t49;
t143 = -t133 * t49 + t73 * t182;
t138 = t48 * t93 - t49 * t94 + t73 * t90 + t74 * t91;
t115 = 0.2e1 * t164;
t109 = t133 * t165;
t108 = t129 * t165;
t105 = -0.2e1 * t156;
t72 = t74 ^ 2;
t61 = -0.2e1 * t107 * t75;
t58 = t133 * t162 + t193;
t57 = t129 * t162 - t189;
t54 = -t146 * qJD(6) - t128 * t108 + t132 * t109;
t53 = -t147 * qJD(6) - t132 * t108 - t128 * t109;
t47 = t106 * t74;
t45 = 0.2e1 * t75 * t106 - 0.2e1 * t107 * t76;
t35 = 0.2e1 * t73 * t49;
t31 = -t106 * t49 - t76 * t73;
t29 = -t74 * t156 + t48 * t185;
t26 = -t148 * qJD(6) - t128 * t57 + t132 * t58;
t25 = t149 * qJD(6) - t128 * t58 - t132 * t57;
t21 = qJD(5) * t157 - t184 * t48;
t14 = t48 * t186 - t128 * t170 - t181 * t194 + (t177 * t190 + t195) * t132;
t13 = -t106 * t48 - t177 * t46;
t9 = t13 * t107 + t47 * t75;
t8 = -t202 * qJD(5) + t161;
t3 = -t13 * t106 - t107 * t14 + t75 * t46 + t47 * t76;
t2 = -t150 * qJD(6) + t166;
t1 = t151 * qJD(6) - t128 * t5 - t209;
t4 = [0, 0, 0, 0.2e1 * t131 * t178, 0.2e1 * (-t131 ^ 2 + t134 ^ 2) * qJD(2), 0, 0, 0, t131 * t176, t134 * t176, -0.2e1 * t77 * t100 - 0.2e1 * t62 * t101 - 0.2e1 * t78 * t140 - 0.2e1 * t63 * t142, 0.2e1 * t168 * t122 + 0.2e1 * t77 * t62 + 0.2e1 * t78 * t63, 0.2e1 * t207, -0.2e1 * t153, 0, 0, 0, 0.2e1 * t87 * t49 + 0.2e1 * t82 * t73, 0.2e1 * t87 * t48 + 0.2e1 * t82 * t74, 0.2e1 * t126 * t207 - 0.2e1 * t72 * t164, 0.2e1 * t156 * t72 + t48 * t157, 0.2e1 * t153 * t133 - 0.2e1 * t73 * t170, -0.2e1 * t153 * t129 - 0.2e1 * t73 * t169, t35, 0.2e1 * t74 * t38 + 0.2e1 * t40 * t49 + 0.2e1 * t8 * t73 + 0.2e1 * (-t44 * t49 + t154) * t129, 0.2e1 * t154 * t133 - 0.2e1 * t202 * t49 - 0.2e1 * t74 * t37 + 0.2e1 * t7 * t73, -0.2e1 * t47 * t13, -0.2e1 * t13 * t46 + 0.2e1 * t47 * t14, 0.2e1 * t13 * t73 - 0.2e1 * t47 * t49, -0.2e1 * t14 * t73 - 0.2e1 * t46 * t49, t35, 0.2e1 * t12 * t46 + 0.2e1 * t32 * t14 - 0.2e1 * t151 * t49 + 0.2e1 * t2 * t73, 0.2e1 * t1 * t73 - 0.2e1 * t12 * t47 + 0.2e1 * t32 * t13 - 0.2e1 * t150 * t49; 0, 0, 0, 0, 0, t178, -t179, 0, -pkin(7) * t178, pkin(7) * t179 (-t188 * t100 - t127 * t140) * pkin(2) (t127 * t63 + t188 * t62) * pkin(2), 0, 0, t48, -t49, 0, -t20, t19, t29, t21, t34, -t143, 0, t37 + (-t152 * qJD(5) - t20) * t133 + t138 * t129, t138 * t133 + t152 * t182 + t203, t9, t3, t30, t31, 0, t86 * t14 - t149 * t49 + t26 * t73 + t81 * t46 + t205, t86 * t13 - t148 * t49 + t25 * t73 - t81 * t47 + t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t91, 0.2e1 * t90, t115, t105, 0, 0, 0, 0.2e1 * t160, 0.2e1 * t199, t61, t45, 0, 0, 0, 0.2e1 * t201, 0.2e1 * t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, 0, 0, 0, 0, 0, t49, t48, 0, 0, 0, 0, 0, -t143, -t34, 0, 0, 0, 0, 0, t31, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t49, 0, -t20, t19, t29, t21, t34, -t143, 0, t37 + (-pkin(4) * t48 - pkin(9) * t49) * t129 + (-t20 + (-pkin(4) * t74 - pkin(9) * t73) * qJD(5)) * t133, pkin(4) * t144 + pkin(9) * t143 + t203, t9, t3, t30, t31, 0, t121 * t14 + t147 * t49 + t173 * t46 + t54 * t73 + t205, t121 * t13 - t146 * t49 - t173 * t47 + t53 * t73 + t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t90, t115, t105, 0, 0, 0, t160 - t175, -t174 + t199, t61, t45, 0, 0, 0, t198 + t201, t197 + t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t105, 0, 0, 0, -0.2e1 * t175, -0.2e1 * t174, t61, t45, 0, 0, 0, 0.2e1 * t198, 0.2e1 * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, -t216, t49, t8, t7, 0, 0, t13, -t14, t49, t132 * t214 + (t167 * t128 - t192) * qJD(6) + t166, -t209 + (-t5 - t214) * t128 + (t132 * t167 + t196) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, -t182, 0, -t94 * t123 + t193, t94 * t182 + t189, 0, 0, -t75, -t76, 0, t26, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, -t123, 0, 0, 0, 0, 0, -t76, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, -t182, 0, -pkin(9) * t123, pkin(9) * t182, 0, 0, -t75, -t76, 0, t54, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t172, -0.2e1 * t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, t49, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t76, 0, t26, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t76, 0, t54, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, -t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t4;
