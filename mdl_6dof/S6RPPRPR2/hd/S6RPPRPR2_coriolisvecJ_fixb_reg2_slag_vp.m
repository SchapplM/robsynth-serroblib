% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:37
% EndTime: 2019-03-09 01:42:43
% DurationCPUTime: 2.10s
% Computational Cost: add. (3930->277), mult. (9294->349), div. (0->0), fcn. (6674->8), ass. (0->156)
t111 = sin(pkin(10));
t113 = cos(pkin(10));
t116 = sin(qJ(4));
t198 = cos(qJ(4));
t96 = t198 * t111 + t116 * t113;
t208 = t96 * qJD(1);
t214 = qJD(6) + t208;
t117 = cos(qJ(6));
t115 = sin(qJ(6));
t163 = t115 * qJD(4);
t159 = t198 * t113;
t148 = qJD(1) * t159;
t172 = t116 * t111;
t158 = qJD(1) * t172;
t86 = -t148 + t158;
t64 = -t117 * t86 + t163;
t156 = t214 * t64;
t164 = qJD(6) * t117;
t91 = t96 * qJD(4);
t76 = qJD(1) * t91;
t42 = qJD(6) * t163 - t115 * t76 - t86 * t164;
t218 = t42 - t156;
t125 = t96 * qJD(3);
t124 = qJD(1) * t125;
t108 = t113 * qJD(2);
t180 = pkin(7) * qJD(1);
t105 = sin(pkin(9)) * pkin(1) + qJ(3);
t98 = t105 * qJD(1);
t62 = t108 + (-t98 - t180) * t111;
t73 = t111 * qJD(2) + t113 * t98;
t63 = t113 * t180 + t73;
t32 = t116 * t62 + t198 * t63;
t18 = t32 * qJD(4) + t124;
t101 = qJD(4) * t148;
t75 = qJD(4) * t158 - t101;
t13 = -t75 * pkin(5) + t18;
t177 = t75 * qJ(5);
t134 = -qJD(5) * t208 + t177;
t201 = pkin(4) + pkin(8);
t16 = t201 * t76 + t134;
t31 = t116 * t63 - t198 * t62;
t169 = -qJD(5) - t31;
t168 = pkin(5) * t208 - t169;
t14 = -t201 * qJD(4) + t168;
t97 = -cos(pkin(9)) * pkin(1) - t113 * pkin(3) - pkin(2);
t82 = qJD(1) * t97 + qJD(3);
t123 = -qJ(5) * t208 + t82;
t24 = t201 * t86 + t123;
t6 = t115 * t14 + t117 * t24;
t2 = -qJD(6) * t6 - t115 * t16 + t117 * t13;
t206 = t214 * t6 + t2;
t135 = t115 * t24 - t117 * t14;
t1 = -qJD(6) * t135 + t115 * t13 + t117 * t16;
t146 = t135 * t214 + t1;
t209 = t117 * t214;
t66 = t117 * qJD(4) + t115 * t86;
t219 = t66 * t209;
t155 = t115 * t214;
t68 = t117 * t75;
t131 = -t155 * t214 - t68;
t95 = -t159 + t172;
t192 = t75 * t95;
t139 = -t208 * t91 + t192;
t157 = qJD(4) * t198;
t165 = qJD(4) * t116;
t90 = t111 * t165 - t113 * t157;
t207 = -t96 * t76 + t90 * t86;
t217 = -t139 + t207;
t216 = t139 + t207;
t215 = t208 * qJD(4);
t202 = t86 ^ 2;
t85 = t208 ^ 2;
t212 = -t202 - t85;
t211 = -t202 + t85;
t54 = t75 * t96;
t210 = -t208 * t90 - t54;
t199 = t86 * pkin(5);
t29 = -qJD(4) * qJ(5) - t32;
t19 = -t29 - t199;
t205 = t19 * t214 + t201 * t75;
t160 = t95 * t164;
t191 = t214 * t91;
t204 = -t115 * (-t191 + t192) + t214 * t160;
t102 = qJD(3) * t159;
t185 = pkin(7) + t105;
t92 = t185 * t111;
t93 = t185 * t113;
t33 = (qJD(3) * t111 + qJD(4) * t93) * t116 + t92 * t157 - t102;
t142 = t42 * t95 - t66 * t91;
t161 = qJD(6) * t115 * t95;
t203 = -t117 * t142 - t66 * t161;
t200 = t76 * pkin(4);
t48 = t116 * t93 + t198 * t92;
t197 = t18 * t48;
t196 = t18 * t95;
t195 = t64 * t91;
t194 = t66 * t64;
t193 = t66 * t86;
t190 = t86 * t64;
t189 = t86 * t208;
t184 = -t42 * t96 - t66 * t90;
t183 = -t29 - t32;
t69 = t117 * t76;
t43 = qJD(6) * t66 - t69;
t178 = t43 * t115;
t176 = t86 * qJ(5);
t174 = qJD(4) * t90;
t173 = qJD(4) * t91;
t171 = t33 * qJD(4);
t49 = -t116 * t92 + t198 * t93;
t34 = qJD(4) * t49 + t125;
t170 = t34 * qJD(4);
t166 = t111 ^ 2 + t113 ^ 2;
t162 = -t115 * t195 - t64 * t160 - t95 * t178;
t152 = -qJD(1) * t102 + qJD(3) * t158 - t62 * t157 + t63 * t165;
t151 = qJD(1) * t166;
t15 = -qJD(4) * qJD(5) + t152;
t11 = -t76 * pkin(5) - t15;
t147 = qJD(6) * t201 * t214 + t11;
t145 = t11 * t95 + t19 * t91;
t144 = t115 * t6 - t117 * t135;
t143 = t115 * t135 + t117 * t6;
t141 = -t96 * t43 + t90 * t64;
t138 = t76 * t95 + t86 * t91;
t137 = t111 * (-t111 * t98 + t108) - t113 * t73;
t127 = -t96 * qJ(5) + t97;
t35 = t201 * t95 + t127;
t38 = t96 * pkin(5) + t48;
t10 = t115 * t38 + t117 * t35;
t9 = -t115 * t35 + t117 * t38;
t133 = t90 * qJ(5) - t96 * qJD(5);
t130 = t117 * t191 - t161 * t214 - t95 * t68;
t128 = t115 * t75 - t209 * t214;
t122 = qJD(6) * t143 + t1 * t115 + t2 * t117;
t121 = t18 * t96 + t208 * t34 + t33 * t86 - t48 * t75 - t49 * t76;
t40 = t86 * pkin(4) + t123;
t120 = t208 * t40 + t124;
t79 = qJD(4) * t86;
t52 = pkin(4) * t208 + t176;
t50 = t75 - t79;
t47 = t95 * pkin(4) + t127;
t41 = t91 * pkin(4) + t133;
t39 = -t95 * pkin(5) + t49;
t37 = t117 * t43;
t36 = t201 * t208 + t176;
t30 = t134 + t200;
t28 = -qJD(4) * pkin(4) - t169;
t25 = t201 * t91 + t133;
t23 = -t90 * pkin(5) + t34;
t22 = -t91 * pkin(5) - t33;
t21 = t32 - t199;
t8 = t115 * t21 + t117 * t36;
t7 = -t115 * t36 + t117 * t21;
t4 = -qJD(6) * t10 - t115 * t25 + t117 * t23;
t3 = qJD(6) * t9 + t115 * t23 + t117 * t25;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t151 (t105 * t151 - t137) * qJD(3), t210, t216, -t174, t138, -t173, 0, t97 * t76 + t82 * t91 - t170, -t97 * t75 - t82 * t90 + t171, t152 * t95 - t31 * t90 - t32 * t91 + t121, -t152 * t49 + t31 * t34 - t32 * t33 + t197, 0, t174, t173, t210, t216, t138, t15 * t95 - t28 * t90 + t29 * t91 + t121, -t30 * t95 - t40 * t91 - t41 * t86 - t47 * t76 + t170, -t208 * t41 - t30 * t96 + t40 * t90 + t47 * t75 - t171, -t15 * t49 + t28 * t34 + t29 * t33 + t30 * t47 + t40 * t41 + t197, -t115 * t142 + t160 * t66, t162 + t203, t184 + t204, t64 * t161 + (-t43 * t95 - t195) * t117, t130 + t141, -t214 * t90 - t54, -t117 * t145 + t135 * t90 + t161 * t19 + t2 * t96 + t214 * t4 + t22 * t64 + t39 * t43 - t9 * t75, -t1 * t96 + t10 * t75 + t115 * t145 + t160 * t19 - t214 * t3 + t22 * t66 - t39 * t42 + t6 * t90, -t10 * t43 - t3 * t64 - t4 * t66 + t9 * t42 + t143 * t91 + (-qJD(6) * t144 + t1 * t117 - t2 * t115) * t95, t1 * t10 + t11 * t39 - t135 * t4 + t19 * t22 + t2 * t9 + t6 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, t174, t217, -t152 * t96 + t31 * t91 - t32 * t90 + t196, 0, 0, 0, 0, 0, 0, t217, t173, -t174, -t15 * t96 + t28 * t91 + t29 * t90 + t196, 0, 0, 0, 0, 0, 0, t130 - t141, t184 - t204, t162 - t203, t11 * t96 + t122 * t95 + t144 * t91 - t19 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166 * qJD(1) ^ 2, t137 * qJD(1), 0, 0, 0, 0, 0, 0, 0.2e1 * t215, t101 + (-t86 - t158) * qJD(4), t212, -t208 * t31 + t32 * t86, 0, 0, 0, 0, 0, 0, t212, -0.2e1 * t215, t75 + t79, t200 + t177 - t29 * t86 + (-qJD(5) - t28) * t208, 0, 0, 0, 0, 0, 0, t128 + t190, t193 - t131, -t218 * t115 + t219 - t37, -t115 * t206 + t146 * t117 + t19 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t211, t101 + (t86 - t158) * qJD(4), -t189, 0, 0 -(qJD(3) + t82) * t208, -t31 * qJD(4) + t82 * t86 + t152, 0, 0, 0, t50, 0, t189, t211, -t189, pkin(4) * t75 - qJ(5) * t76 + t183 * t208 + (t28 + t169) * t86, t52 * t86 + t120, -t40 * t86 + t52 * t208 + (0.2e1 * qJD(5) + t31) * qJD(4) - t152, -t18 * pkin(4) - t15 * qJ(5) + t169 * t29 - t28 * t32 - t40 * t52, -t117 * t42 - t155 * t66, -t37 - t219 + (t42 + t156) * t115, t131 + t193, t209 * t64 + t178, t128 - t190, t214 * t86, qJ(5) * t43 + t147 * t115 + t117 * t205 - t135 * t86 + t168 * t64 - t214 * t7, -qJ(5) * t42 - t115 * t205 + t147 * t117 + t168 * t66 + t214 * t8 - t6 * t86, t8 * t64 + t7 * t66 + (-t201 * t42 - t6 * t208 - t2 + (t201 * t64 - t6) * qJD(6)) * t117 + (t201 * t43 - t135 * t208 - t1 + (-t201 * t66 - t135) * qJD(6)) * t115, t11 * qJ(5) - t122 * t201 + t135 * t7 + t168 * t19 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t189, -qJD(4) ^ 2 - t85, -qJD(4) * t183 + t120, 0, 0, 0, 0, 0, 0, -qJD(4) * t64 + t131, -qJD(4) * t66 + t128, t218 * t117 + (t214 * t66 - t43) * t115, -t19 * qJD(4) + t146 * t115 + t117 * t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, -t64 ^ 2 + t66 ^ 2, -t218, -t194, t69 + (-qJD(6) + t214) * t66, -t75, -t19 * t66 + t206, t19 * t64 - t146, 0, 0;];
tauc_reg  = t5;
