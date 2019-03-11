% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:10
% EndTime: 2019-03-09 01:56:16
% DurationCPUTime: 1.96s
% Computational Cost: add. (3972->298), mult. (8608->372), div. (0->0), fcn. (5952->6), ass. (0->165)
t104 = cos(qJ(6));
t101 = cos(pkin(9));
t105 = cos(qJ(4));
t166 = t105 * t101;
t152 = qJD(1) * t166;
t100 = sin(pkin(9));
t194 = sin(qJ(4));
t153 = t194 * t100;
t89 = qJD(1) * t153;
t73 = -t89 + t152;
t66 = qJD(6) + t73;
t141 = t104 * t66;
t103 = sin(qJ(6));
t78 = t105 * t100 + t194 * t101;
t169 = qJD(1) * t78;
t50 = t104 * qJD(4) + t103 * t169;
t214 = t50 * t141;
t102 = -pkin(1) - qJ(3);
t207 = t102 * qJD(1);
t87 = qJD(2) + t207;
t148 = -pkin(7) * qJD(1) + t87;
t61 = t148 * t100;
t62 = t148 * t101;
t40 = -t105 * t62 + t194 * t61;
t162 = -qJD(5) - t40;
t161 = t73 * pkin(5) - t162;
t199 = pkin(4) + pkin(8);
t16 = -t199 * qJD(4) + t161;
t99 = qJD(1) * qJ(2);
t94 = qJD(3) + t99;
t95 = t100 * pkin(3);
t82 = qJD(1) * t95 + t94;
t119 = -t73 * qJ(5) + t82;
t23 = t169 * t199 + t119;
t121 = t103 * t23 - t104 * t16;
t151 = qJD(3) * t194;
t136 = qJD(1) * t151;
t160 = qJD(3) * t105;
t149 = qJD(1) * t160;
t117 = -t100 * t136 + t101 * t149;
t41 = t105 * t61 + t194 * t62;
t21 = t41 * qJD(4) + t117;
t75 = t78 * qJD(4);
t59 = qJD(1) * t75;
t13 = -t59 * pkin(5) + t21;
t98 = qJD(1) * qJD(2);
t147 = t59 * qJ(5) + t98;
t115 = -t73 * qJD(5) + t147;
t86 = qJD(4) * t89;
t60 = qJD(4) * t152 - t86;
t14 = t199 * t60 + t115;
t1 = -t121 * qJD(6) + t103 * t13 + t104 * t14;
t134 = t121 * t66 + t1;
t6 = t103 * t16 + t104 * t23;
t2 = -qJD(6) * t6 - t103 * t14 + t104 * t13;
t213 = t6 * t66 + t2;
t200 = t169 ^ 2;
t69 = t73 ^ 2;
t212 = -t200 - t69;
t211 = -t200 + t69;
t157 = t103 * qJD(4);
t48 = -t104 * t169 + t157;
t145 = t48 * t66;
t158 = qJD(6) * t104;
t31 = qJD(6) * t157 - t103 * t60 - t158 * t169;
t210 = t31 - t145;
t32 = t50 * qJD(6) - t104 * t60;
t209 = t50 * t66 - t32;
t79 = -t153 + t166;
t45 = t79 * t59;
t124 = -t75 * t73 - t45;
t167 = qJD(4) * t75;
t208 = qJD(1) * t169 + t167;
t176 = t100 ^ 2 + t101 ^ 2;
t206 = t176 * qJD(3);
t150 = qJD(4) * t194;
t159 = qJD(4) * t105;
t144 = t100 * t149 + t101 * t136 + t61 * t150 - t62 * t159;
t17 = -qJD(4) * qJD(5) + t144;
t33 = -qJD(4) * pkin(4) - t162;
t36 = -qJD(4) * qJ(5) - t41;
t76 = -t100 * t150 + t101 * t159;
t205 = t17 * t78 - t33 * t75 + t36 * t76;
t204 = t144 * t78 - t40 * t75 - t41 * t76;
t195 = t169 * pkin(5);
t22 = -t36 - t195;
t203 = t199 * t59 + t22 * t66;
t178 = -pkin(7) + t102;
t80 = t178 * t100;
t81 = t178 * t101;
t34 = t78 * qJD(3) + t80 * t150 - t81 * t159;
t202 = qJD(4) * (t73 + t152) - t86;
t201 = qJD(4) * (-t73 + t152) - t86;
t156 = 0.2e1 * t98;
t196 = t60 * pkin(4);
t46 = -t105 * t81 + t194 * t80;
t191 = t21 * t46;
t190 = t21 * t79;
t185 = t50 * t48;
t183 = t50 * t169;
t182 = t59 * t78;
t181 = t169 * t48;
t180 = t169 * t73;
t177 = -t36 - t41;
t175 = t103 * t59;
t174 = t104 * t31;
t54 = t104 * t59;
t172 = t32 * t103;
t171 = t169 * qJ(5);
t91 = qJ(2) + t95;
t168 = qJD(4) * t169;
t165 = t34 * qJD(4);
t47 = t105 * t80 + t194 * t81;
t35 = qJD(4) * t47 - t100 * t151 + t101 * t160;
t164 = t35 * qJD(4);
t163 = t76 * qJD(4);
t155 = qJD(6) * t103 * t78;
t154 = t78 * t158;
t146 = qJD(1) * t176;
t143 = t103 * t66;
t140 = qJD(6) * t79 + qJD(1);
t139 = -t1 * t79 + t6 * t75;
t138 = -t121 * t75 - t2 * t79;
t137 = -t79 * qJ(5) + t91;
t11 = -t60 * pkin(5) - t17;
t135 = qJD(6) * t199 * t66 + t11;
t131 = t103 * t121 + t104 * t6;
t130 = t11 * t78 + t22 * t76;
t129 = -t31 * t79 - t50 * t75;
t128 = -t78 * t31 + t76 * t50;
t127 = t32 * t79 - t48 * t75;
t126 = -t78 * t32 - t76 * t48;
t125 = t66 * t76 - t182;
t123 = t169 * t76 + t60 * t78;
t30 = t199 * t78 + t137;
t38 = t79 * pkin(5) + t46;
t10 = t103 * t38 + t104 * t30;
t9 = -t103 * t30 + t104 * t38;
t120 = qJD(1) * t73 + t163;
t116 = -t66 * t143 - t54;
t114 = t75 * qJ(5) - t79 * qJD(5) + qJD(2);
t112 = -t141 * t66 + t175;
t37 = pkin(4) * t169 + t119;
t111 = t37 * t73 + t117;
t110 = -t123 - t124;
t109 = t169 * t75 - t79 * t60 - t73 * t76 + t182;
t108 = t169 * t34 + t35 * t73 - t46 * t59 - t47 * t60 + t190;
t107 = qJD(1) ^ 2;
t44 = t73 * pkin(4) + t171;
t43 = t78 * pkin(4) + t137;
t42 = -t168 + t59;
t39 = -t78 * pkin(5) + t47;
t29 = t104 * t32;
t28 = t76 * pkin(4) + t114;
t27 = t199 * t73 + t171;
t26 = t41 - t195;
t24 = t115 + t196;
t19 = -t75 * pkin(5) + t35;
t18 = -t76 * pkin(5) - t34;
t15 = t199 * t76 + t114;
t8 = t103 * t26 + t104 * t27;
t7 = -t103 * t27 + t104 * t26;
t4 = -t10 * qJD(6) - t103 * t15 + t104 * t19;
t3 = t9 * qJD(6) + t103 * t19 + t104 * t15;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, qJ(2) * t156, 0, 0, 0, 0, 0, 0, t100 * t156, t101 * t156, 0.2e1 * qJD(3) * t146 (t94 + t99) * qJD(2) + (-t87 - t207) * t206, t124, t109, -t167, t123, -t163, 0, 0.2e1 * t169 * qJD(2) + t91 * t60 + t82 * t76 - t164, t165 - t91 * t59 - t82 * t75 + (qJD(1) * t79 + t73) * qJD(2), t108 + t204, -t144 * t47 + t191 - t41 * t34 + t40 * t35 + (qJD(1) * t91 + t82) * qJD(2), 0, t167, t163, t124, t109, t123, t108 + t205, -t169 * t28 - t24 * t78 - t37 * t76 - t43 * t60 + t164, -t24 * t79 - t28 * t73 + t37 * t75 + t43 * t59 - t165, -t17 * t47 + t24 * t43 + t37 * t28 + t33 * t35 + t36 * t34 + t191, t103 * t128 + t154 * t50 (-t103 * t48 + t104 * t50) * t76 + (-t172 - t174 + (-t103 * t50 - t104 * t48) * qJD(6)) * t78, t103 * t125 + t154 * t66 + t129, t104 * t126 + t155 * t48, t104 * t125 - t155 * t66 - t127, -t66 * t75 - t45, -t104 * t130 + t155 * t22 + t18 * t48 + t39 * t32 + t4 * t66 - t9 * t59 - t138, t10 * t59 + t103 * t130 + t154 * t22 + t18 * t50 - t3 * t66 - t39 * t31 + t139, -t10 * t32 - t3 * t48 + t9 * t31 - t4 * t50 + t131 * t76 + (t1 * t104 - t2 * t103 + (-t103 * t6 + t104 * t121) * qJD(6)) * t78, t1 * t10 + t11 * t39 - t121 * t4 + t22 * t18 + t2 * t9 + t6 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t107 * qJ(2), 0, 0, 0, 0, 0, 0, -t107 * t100, -t107 * t101, 0 (-t94 - t206) * qJD(1), 0, 0, 0, 0, 0, 0, -t208, -t120, t110, -t82 * qJD(1) - t190 - t204, 0, 0, 0, 0, 0, 0, t110, t208, t120, -t37 * qJD(1) - t190 - t205, 0, 0, 0, 0, 0, 0, t79 * t54 + (t103 * t140 + t104 * t75) * t66 - t126, -t79 * t175 + (-t103 * t75 + t104 * t140) * t66 + t128 (t140 * t48 + t129) * t104 + (-t140 * t50 + t127) * t103 (-t140 * t6 + t138) * t104 + (-t121 * t140 + t139) * t103 + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176 * t107, t146 * t87 + t98, 0, 0, 0, 0, 0, 0, t202, -0.2e1 * t168, t212, t169 * t41 - t40 * t73 + t98, 0, 0, 0, 0, 0, 0, t212, -t202, t168 + t59, t196 - t36 * t169 + (-qJD(5) - t33) * t73 + t147, 0, 0, 0, 0, 0, 0, t112 + t181, t66 ^ 2 * t103 + t183 + t54, -t210 * t103 + t214 - t29, -t213 * t103 + t134 * t104 + t169 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t211, 0, -t180, -t201, 0, -t82 * t73 - t117, -t40 * qJD(4) + t169 * t82 + t144, 0, 0, 0, t42, t201, t180, t211, -t180, pkin(4) * t59 - qJ(5) * t60 + t177 * t73 + (t33 + t162) * t169, t169 * t44 + t111, -t37 * t169 + t44 * t73 + (0.2e1 * qJD(5) + t40) * qJD(4) - t144, -t21 * pkin(4) - t17 * qJ(5) + t162 * t36 - t33 * t41 - t37 * t44, -t143 * t50 - t174, -t29 - t214 + (t31 + t145) * t103, t116 + t183, t141 * t48 + t172, t112 - t181, t66 * t169, qJ(5) * t32 + t135 * t103 + t203 * t104 - t121 * t169 + t161 * t48 - t7 * t66, -qJ(5) * t31 - t203 * t103 + t135 * t104 + t161 * t50 - t169 * t6 + t8 * t66, t8 * t48 + t7 * t50 + (-t199 * t31 - t6 * t73 - t2 + (t199 * t48 - t6) * qJD(6)) * t104 + (t199 * t32 - t121 * t73 - t1 + (-t199 * t50 - t121) * qJD(6)) * t103, t11 * qJ(5) + t121 * t7 - t6 * t8 + t161 * t22 - (qJD(6) * t131 + t1 * t103 + t2 * t104) * t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t180, -qJD(4) ^ 2 - t69, -t177 * qJD(4) + t111, 0, 0, 0, 0, 0, 0, -qJD(4) * t48 + t116, -qJD(4) * t50 + t112, t209 * t103 + t210 * t104, -t22 * qJD(4) + t134 * t103 + t213 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, -t48 ^ 2 + t50 ^ 2, -t210, -t185, t209, -t59, -t22 * t50 + t213, t22 * t48 - t134, 0, 0;];
tauc_reg  = t5;
