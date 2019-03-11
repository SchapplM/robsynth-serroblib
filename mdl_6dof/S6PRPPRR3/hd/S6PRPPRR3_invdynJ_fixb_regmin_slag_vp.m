% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPPRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:06
% EndTime: 2019-03-08 19:23:11
% DurationCPUTime: 2.12s
% Computational Cost: add. (1562->304), mult. (3258->429), div. (0->0), fcn. (2670->12), ass. (0->167)
t119 = cos(qJ(5));
t96 = t119 * qJD(2) + qJD(6);
t221 = t96 - qJD(6);
t108 = sin(pkin(11));
t111 = cos(pkin(11));
t109 = sin(pkin(10));
t112 = cos(pkin(10));
t117 = sin(qJ(2));
t113 = cos(pkin(6));
t120 = cos(qJ(2));
t193 = t113 * t120;
t65 = t109 * t193 + t112 * t117;
t194 = t113 * t117;
t66 = -t109 * t194 + t112 * t120;
t155 = t66 * t108 - t65 * t111;
t63 = t109 * t117 - t112 * t193;
t64 = t109 * t120 + t112 * t194;
t156 = t64 * t108 - t63 * t111;
t110 = sin(pkin(6));
t56 = (t108 * t117 + t111 * t120) * t110;
t135 = g(1) * t155 + g(2) * t156 + g(3) * t56;
t187 = qJD(1) * t110;
t165 = t120 * t187;
t166 = t117 * t187;
t46 = t108 * t165 - t111 * t166;
t154 = t108 * qJD(3) - t46;
t220 = t154 * qJD(2) - t135;
t122 = qJD(5) ^ 2;
t121 = -pkin(2) - pkin(3);
t87 = qJD(2) * t166;
t176 = qJDD(1) * t110;
t90 = t120 * t176;
t167 = qJDD(3) + t87 - t90;
t43 = t121 * qJDD(2) + t167;
t175 = qJDD(2) * qJ(3);
t89 = t117 * t176;
t44 = t175 + t89 + (qJD(3) + t165) * qJD(2);
t9 = -t108 * t44 + t111 * t43;
t7 = qJDD(2) * pkin(4) - t9;
t78 = -t108 * qJ(3) + t111 * t121;
t74 = pkin(4) - t78;
t79 = t111 * qJ(3) + t108 * t121;
t75 = -pkin(8) + t79;
t219 = qJDD(2) * t74 - t122 * t75 + t220 + t7;
t103 = t119 * qJDD(2);
t116 = sin(qJ(5));
t178 = qJD(2) * qJD(5);
t161 = t116 * t178;
t218 = -t161 + t103;
t196 = t110 * t120;
t198 = t110 * t117;
t217 = t108 * t196 - t111 * t198;
t197 = t110 * t119;
t27 = t63 * t108 + t64 * t111;
t31 = t65 * t108 + t66 * t111;
t38 = t113 * t119 - t116 * t217;
t137 = g(1) * (-t109 * t197 - t31 * t116) + g(2) * (t112 * t197 - t27 * t116) - g(3) * t38;
t148 = -pkin(5) * t116 + pkin(9) * t119;
t91 = -t113 * qJDD(1) + qJDD(4);
t205 = t119 * t91;
t147 = qJD(3) - t165;
t67 = t121 * qJD(2) + t147;
t77 = qJD(2) * qJ(3) + t166;
t37 = t108 * t67 + t111 * t77;
t33 = -qJD(2) * pkin(8) + t37;
t93 = -t113 * qJD(1) + qJD(4);
t22 = t116 * t93 + t119 * t33;
t10 = t108 * t43 + t111 * t44;
t8 = -qJDD(2) * pkin(8) + t10;
t2 = -qJDD(5) * pkin(5) + t22 * qJD(5) + t116 * t8 - t205;
t216 = (pkin(9) * qJD(6) + t148 * qJD(2)) * t96 + t137 + t2;
t149 = t119 * pkin(5) + t116 * pkin(9);
t50 = t149 + t74;
t70 = -qJDD(6) - t218;
t215 = (-qJD(5) * t148 - t154) * t96 + t50 * t70;
t16 = qJD(5) * pkin(9) + t22;
t210 = t75 * t96;
t214 = -qJD(6) * (t16 + t210) - t135;
t115 = sin(qJ(6));
t186 = qJD(2) * t116;
t164 = t115 * t186;
t118 = cos(qJ(6));
t179 = t118 * qJD(5);
t72 = t164 + t179;
t212 = t72 * t96;
t180 = t115 * qJD(5);
t73 = t118 * t186 - t180;
t211 = t73 * t96;
t208 = pkin(2) * t196 + qJ(3) * t198;
t207 = t116 * t91;
t206 = t118 * t96;
t162 = t119 * t178;
t172 = t116 * qJDD(2);
t134 = t162 + t172;
t177 = qJD(5) * qJD(6);
t34 = qJD(6) * t164 + t115 * qJDD(5) + (-t134 + t177) * t118;
t204 = t34 * t115;
t203 = qJD(5) * t72;
t202 = qJD(5) * t73;
t201 = qJD(5) * t75;
t200 = qJDD(2) * pkin(2);
t199 = t110 * t116;
t195 = t111 * t116;
t192 = t115 * t119;
t191 = t118 * t119;
t190 = qJDD(1) - g(3);
t105 = t116 ^ 2;
t189 = -t119 ^ 2 + t105;
t123 = qJD(2) ^ 2;
t188 = t122 + t123;
t185 = qJD(5) * t116;
t184 = qJD(5) * t119;
t183 = qJD(6) * t115;
t182 = qJD(6) * t118;
t174 = qJDD(5) * t116;
t173 = qJDD(5) * t119;
t169 = t96 * t180;
t168 = t96 * t179;
t159 = -t63 * pkin(2) + t64 * qJ(3);
t158 = -t65 * pkin(2) + t66 * qJ(3);
t36 = -t108 * t77 + t111 * t67;
t32 = qJD(2) * pkin(4) - t36;
t157 = qJD(2) * t32 - t8;
t49 = qJD(1) * t56;
t153 = t111 * qJD(3) - t49;
t152 = 0.2e1 * t161;
t35 = t118 * qJDD(5) + (t116 * t182 + t119 * t180) * qJD(2) + (t172 - t177) * t115;
t151 = -t35 + t169;
t150 = t34 + t168;
t146 = t36 * t108 - t37 * t111;
t23 = qJD(2) * t149 + t32;
t4 = t115 * t23 + t118 * t16;
t145 = t115 * t16 - t118 * t23;
t39 = -t113 * t116 - t119 * t217;
t14 = t56 * t115 + t39 * t118;
t13 = -t39 * t115 + t56 * t118;
t144 = t116 * t33 - t119 * t93;
t143 = t117 * (-qJD(2) * pkin(2) + t147) + t120 * t77;
t142 = -g(1) * t65 - g(2) * t63 + g(3) * t196;
t140 = -t115 * t70 + t182 * t96;
t139 = t118 * t70 + t183 * t96;
t138 = -t142 + t90;
t136 = -g(1) * t31 - g(2) * t27 + g(3) * t217;
t132 = -qJDD(3) + t138;
t131 = g(1) * t66 + g(2) * t64 + g(3) * t198 - t89;
t130 = -t140 - t203;
t129 = t139 - t202;
t128 = -qJD(6) * t50 * t96 - t136;
t15 = -qJD(5) * pkin(5) + t144;
t127 = pkin(9) * t70 + (t15 - t144) * t96;
t126 = -qJDD(5) * t75 + (-qJD(2) * t74 - t153 - t32) * qJD(5);
t1 = qJDD(5) * pkin(9) - qJD(5) * t144 + t119 * t8 + t207;
t124 = -t15 * qJD(5) - qJD(6) * t23 - t153 * t96 + t75 * t70 - t1;
t83 = -t122 * t116 + t173;
t82 = -t122 * t119 - t174;
t62 = (qJDD(2) * t117 + t120 * t123) * t110;
t61 = (qJDD(2) * t120 - t117 * t123) * t110;
t48 = qJD(2) * t56;
t47 = t217 * qJD(2);
t45 = t167 - t200;
t20 = -t109 * t199 + t31 * t119;
t18 = t112 * t199 + t27 * t119;
t12 = -qJD(5) * t38 + t48 * t119;
t11 = qJD(5) * t39 + t48 * t116;
t6 = t218 * pkin(5) + t134 * pkin(9) + t7;
t5 = t118 * t6;
t3 = [t190, 0, t61, -t62, t61, t62, t113 ^ 2 * qJDD(1) - g(3) + (qJD(2) * t143 + t117 * t44 - t120 * t45) * t110, qJD(2) * t47 + qJDD(2) * t56, qJD(2) * t48 - qJDD(2) * t217, -t10 * t217 - t91 * t113 - t36 * t47 + t37 * t48 - t9 * t56 - g(3), 0, 0, 0, 0, 0, t56 * t103 - t11 * qJD(5) - t38 * qJDD(5) + (t119 * t47 - t185 * t56) * qJD(2), -t56 * t172 - t12 * qJD(5) - t39 * qJDD(5) + (-t116 * t47 - t184 * t56) * qJD(2), 0, 0, 0, 0, 0 (-qJD(6) * t14 - t12 * t115 + t47 * t118) * t96 - t13 * t70 - t11 * t72 - t38 * t35 -(qJD(6) * t13 + t47 * t115 + t12 * t118) * t96 + t14 * t70 - t11 * t73 + t38 * t34; 0, qJDD(2), t138, t131, t132 + 0.2e1 * t200, 0.2e1 * qJD(2) * qJD(3) - t131 + 0.2e1 * t175, -t45 * pkin(2) - g(1) * t158 - g(2) * t159 - g(3) * t208 + t44 * qJ(3) + t77 * qJD(3) - t143 * t187, -t78 * qJDD(2) + t220 - t9, qJD(2) * t153 + t79 * qJDD(2) + t10 + t136, t10 * t79 + t9 * t78 - t37 * t49 + t36 * t46 - g(1) * (-t65 * pkin(3) + t158) - g(2) * (-t63 * pkin(3) + t159) - g(3) * (pkin(3) * t196 + t208) - t146 * qJD(3), t105 * qJDD(2) + t119 * t152, 0.2e1 * t116 * t103 - 0.2e1 * t189 * t178, t82, -t83, 0, t126 * t116 + t219 * t119, -t219 * t116 + t126 * t119, t73 * t119 * t179 + (-t34 * t118 - t183 * t73) * t116 (-t115 * t73 - t118 * t72) * t184 + (t204 - t118 * t35 + (t115 * t72 - t118 * t73) * qJD(6)) * t116 (t34 - t168) * t119 + (t139 + t202) * t116 (t35 + t169) * t119 + (t140 - t203) * t116, -t70 * t119 - t185 * t96, -t215 * t118 + t128 * t115 + (t124 * t115 + t214 * t118 - t72 * t201 + t5) * t119 + (-t15 * t182 - t2 * t115 - t75 * t35 - t153 * t72 + (t115 * t210 + t145) * qJD(5)) * t116, t215 * t115 + t128 * t118 + (-t73 * t201 + t124 * t118 + (-t214 - t6) * t115) * t119 + (t15 * t183 - t2 * t118 + t75 * t34 - t153 * t73 + (t75 * t206 + t4) * qJD(5)) * t116; 0, 0, 0, 0, -qJDD(2), -t123, -t77 * qJD(2) - t132 - t200 + t87, -t111 * qJDD(2) - t108 * t123, t108 * qJDD(2) - t111 * t123, qJD(2) * t146 + t10 * t108 + t9 * t111 + t142, 0, 0, 0, 0, 0 (t152 - t103) * t111 + (-t188 * t119 - t174) * t108 (0.2e1 * t162 + t172) * t111 + (t188 * t116 - t173) * t108, 0, 0, 0, 0, 0, t139 * t111 + (t116 * t151 + t119 * t130) * t108 + (-(t108 * t118 - t111 * t192) * t96 + t72 * t195) * qJD(2), t140 * t111 + (t116 * t150 + t119 * t129) * t108 + ((t108 * t115 + t111 * t191) * t96 + t73 * t195) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) - t190 * t113 + (g(1) * t109 - g(2) * t112) * t110, 0, 0, 0, 0, 0, t83, t82, 0, 0, 0, 0, 0, t116 * t130 - t119 * t151, t116 * t129 - t119 * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 * t123 * t119, t189 * t123, -t172, -t103, qJDD(5), t116 * t157 - t137 + t205, g(1) * t20 + g(2) * t18 + g(3) * t39 + t157 * t119 - t207, -t73 * t206 + t204 (t34 + t212) * t118 + (t35 + t211) * t115 (-t116 * t73 + t96 * t191) * qJD(2) + t140 (t116 * t72 - t96 * t192) * qJD(2) - t139, t96 * t186, pkin(5) * t35 + t127 * t115 - t216 * t118 - t145 * t186 + t22 * t72, -pkin(5) * t34 + t216 * t115 + t127 * t118 - t4 * t186 + t22 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t72, -t72 ^ 2 + t73 ^ 2, t34 - t212, t35 - t211, -t70, -t115 * t1 + t5 + t15 * t73 - g(1) * (-t20 * t115 + t118 * t155) - g(2) * (-t18 * t115 + t118 * t156) - g(3) * t13 + t221 * t4, -t118 * t1 - t115 * t6 - t15 * t72 - g(1) * (-t115 * t155 - t20 * t118) - g(2) * (-t115 * t156 - t18 * t118) + g(3) * t14 - t221 * t145;];
tau_reg  = t3;
