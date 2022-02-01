% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:59
% EndTime: 2022-01-23 09:17:03
% DurationCPUTime: 1.91s
% Computational Cost: add. (2005->271), mult. (5083->385), div. (0->0), fcn. (4002->14), ass. (0->166)
t134 = cos(pkin(8));
t187 = t134 * qJD(1);
t227 = qJD(4) - t187;
t109 = -qJD(5) - t227;
t135 = sin(qJ(5));
t138 = cos(qJ(5));
t132 = sin(pkin(8));
t131 = sin(pkin(9));
t133 = cos(pkin(9));
t136 = sin(qJ(4));
t139 = cos(qJ(4));
t94 = t131 * t139 + t133 * t136;
t149 = qJD(1) * t94;
t66 = t132 * t149;
t195 = qJD(1) * t132;
t178 = t131 * t195;
t164 = t136 * t178;
t201 = t139 * t133;
t181 = t132 * t201;
t69 = qJD(1) * t181 - t164;
t23 = t135 * t69 + t138 * t66;
t214 = t109 * t23;
t189 = qJD(5) * t138;
t190 = qJD(5) * t135;
t184 = qJDD(1) * t131;
t89 = t94 * qJD(4);
t34 = qJDD(1) * t181 + (-qJD(1) * t89 - t136 * t184) * t132;
t191 = qJD(4) * t139;
t177 = t133 * t191;
t35 = -qJD(4) * t164 + (qJD(1) * t177 + qJDD(1) * t94) * t132;
t4 = -t135 * t35 + t138 * t34 - t189 * t66 - t190 * t69;
t230 = t4 - t214;
t156 = -t135 * t66 + t138 * t69;
t229 = t156 * t23;
t137 = sin(qJ(1));
t140 = cos(qJ(1));
t222 = g(1) * t137 - g(2) * t140;
t212 = t156 * t109;
t5 = qJD(5) * t156 + t135 * t34 + t138 * t35;
t228 = -t5 - t212;
t226 = t156 ^ 2 - t23 ^ 2;
t148 = -t133 * t132 * pkin(6) + (-qJ(2) * t131 - pkin(3)) * t134;
t97 = pkin(2) * t134 + qJ(3) * t132 + pkin(1);
t84 = -qJD(1) * t97 + qJD(2);
t74 = t133 * t84;
t36 = qJD(1) * t148 + t74;
t179 = qJ(2) * t187;
t48 = t131 * t84 + t133 * t179;
t41 = -pkin(6) * t178 + t48;
t155 = -t136 * t36 - t139 * t41;
t11 = -pkin(7) * t66 - t155;
t130 = pkin(9) + qJ(4);
t124 = qJ(5) + t130;
t119 = cos(t124);
t182 = t134 * qJDD(1);
t112 = -qJDD(4) + t182;
t185 = qJD(1) * qJD(2);
t175 = t134 * t185;
t193 = qJD(3) * t132;
t57 = -qJD(1) * t193 - qJDD(1) * t97 + qJDD(2);
t52 = t133 * t57;
t19 = qJDD(1) * t148 - t131 * t175 + t52;
t183 = qJDD(1) * t132;
t174 = t131 * t183;
t173 = t133 * t182;
t33 = qJ(2) * t173 + t131 * t57 + t133 * t175;
t21 = -pkin(6) * t174 + t33;
t170 = -t136 * t21 + t139 * t19;
t147 = qJD(4) * t155 + t170;
t2 = -t112 * pkin(4) - t34 * pkin(7) + t147;
t218 = g(3) * t132;
t110 = qJ(2) * t195 + qJD(3);
t85 = pkin(3) * t178 + t110;
t39 = pkin(4) * t66 + t85;
t118 = sin(t124);
t200 = t140 * t118;
t204 = t137 * t119;
t62 = -t134 * t204 + t200;
t199 = t140 * t119;
t205 = t137 * t118;
t64 = t134 * t199 + t205;
t9 = t11 * t190;
t225 = t39 * t23 + t119 * t218 + t9 + g(1) * t64 + (t109 * t11 - t2) * t135 - g(2) * t62;
t209 = qJDD(1) * pkin(1);
t223 = t209 + t222;
t186 = qJ(2) * qJDD(1);
t213 = t11 * t138;
t169 = -t136 * t41 + t139 * t36;
t10 = -pkin(7) * t69 + t169;
t8 = pkin(4) * t227 + t10;
t158 = -t135 * t8 - t213;
t192 = qJD(4) * t136;
t152 = -t136 * t19 - t139 * t21 - t191 * t36 + t192 * t41;
t3 = -pkin(7) * t35 - t152;
t180 = -t135 * t3 + t138 * t2;
t61 = t134 * t205 + t199;
t63 = -t134 * t200 + t204;
t221 = -g(1) * t63 + g(2) * t61 + qJD(5) * t158 + t118 * t218 - t156 * t39 + t180;
t92 = t133 * t97;
t44 = -t92 + t148;
t208 = t131 * t132;
t60 = qJ(2) * t133 * t134 - t131 * t97;
t49 = -pkin(6) * t208 + t60;
t217 = t136 * t44 + t139 * t49;
t216 = -t134 * t149 + t89;
t93 = -t131 * t136 + t201;
t215 = t227 * t93;
t211 = t66 * t227;
t210 = t69 * t227;
t207 = t131 * t134;
t141 = qJD(1) ^ 2;
t206 = t134 * t141;
t121 = sin(t130);
t203 = t137 * t121;
t122 = cos(t130);
t202 = t137 * t122;
t198 = t140 * t121;
t197 = t140 * t122;
t95 = pkin(3) * t208 + qJ(2) * t132;
t128 = t132 ^ 2;
t196 = t134 ^ 2 + t128;
t194 = qJD(2) * t134;
t188 = t132 * qJD(2);
t176 = qJD(5) * t8 + t3;
t90 = qJ(2) * t183 + t132 * t185 + qJDD(3);
t168 = -t136 * t49 + t139 * t44;
t86 = -t131 * t194 - t133 * t193;
t87 = -t131 * t193 + t133 * t194;
t167 = -t136 * t87 + t139 * t86;
t166 = t196 * t141;
t165 = 0.2e1 * t196;
t58 = pkin(3) * t174 + t90;
t162 = qJD(5) * t94 + t216;
t161 = qJD(5) * t93 + t215;
t160 = g(1) * t140 + g(2) * t137;
t82 = t94 * t132;
t83 = t93 * t132;
t37 = t135 * t83 + t138 * t82;
t38 = -t135 * t82 + t138 * t83;
t153 = t185 + t186;
t151 = t136 * t86 + t139 * t87 + t191 * t44 - t192 * t49;
t150 = t165 * t185;
t144 = t128 * t153 + t90 * t132 - t160;
t126 = t140 * qJ(2);
t125 = t137 * qJ(2);
t120 = qJDD(2) - t209;
t106 = -qJDD(5) + t112;
t78 = t134 * t197 + t203;
t77 = -t134 * t198 + t202;
t76 = -t134 * t202 + t198;
t75 = t134 * t203 + t197;
t72 = t132 * t177 - t192 * t208;
t71 = t132 * t89;
t59 = -qJ(2) * t207 - t92;
t50 = pkin(4) * t72 + t188;
t47 = -t131 * t179 + t74;
t46 = pkin(4) * t82 + t95;
t32 = -t153 * t207 + t52;
t16 = pkin(4) * t35 + t58;
t15 = -pkin(7) * t82 + t217;
t14 = -pkin(4) * t134 - pkin(7) * t83 + t168;
t13 = qJD(5) * t38 - t135 * t71 + t138 * t72;
t12 = -qJD(5) * t37 - t135 * t72 - t138 * t71;
t7 = t71 * pkin(7) - qJD(4) * t217 + t167;
t6 = -pkin(7) * t72 + t151;
t1 = [qJDD(1), t222, t160, (-t120 + t223) * t134, t165 * t186 + t150 - t160, -t120 * pkin(1) - g(1) * (-pkin(1) * t137 + t126) - g(2) * (pkin(1) * t140 + t125) + (t186 * t196 + t150) * qJ(2), (-t86 * qJD(1) - t59 * qJDD(1) + t133 * t222 - t32) * t134 + t144 * t131, (t87 * qJD(1) + t60 * qJDD(1) - t131 * t222 + t33) * t134 + t144 * t133, t33 * t60 + t48 * t87 + t32 * t59 + t47 * t86 - g(1) * (-t137 * t97 + t126) - g(2) * (t140 * t97 + t125) + (qJ(2) * t90 + qJD(2) * t110) * t132, t34 * t83 - t69 * t71, -t34 * t82 - t35 * t83 + t66 * t71 - t69 * t72, -t112 * t83 - t134 * t34 - t227 * t71, t112 * t82 + t134 * t35 - t227 * t72, t112 * t134, t167 * t227 - t168 * t112 - t170 * t134 + t66 * t188 + t95 * t35 + t58 * t82 + t85 * t72 - g(1) * t76 - g(2) * t78 + (-t134 * t155 - t217 * t227) * qJD(4), -g(1) * t75 - g(2) * t77 + t112 * t217 - t134 * t152 - t151 * t227 + t188 * t69 + t95 * t34 + t58 * t83 - t85 * t71, t12 * t156 + t38 * t4, -t12 * t23 - t13 * t156 - t37 * t4 - t38 * t5, -t106 * t38 - t109 * t12 - t134 * t4, t106 * t37 + t109 * t13 + t134 * t5, t106 * t134, -(-t135 * t6 + t138 * t7) * t109 - (-t135 * t15 + t138 * t14) * t106 - t180 * t134 + t50 * t23 + t46 * t5 + t16 * t37 + t39 * t13 - g(1) * t62 - g(2) * t64 + (-(-t135 * t14 - t138 * t15) * t109 - t158 * t134) * qJD(5), -g(1) * t61 - g(2) * t63 + t39 * t12 - t9 * t134 + t16 * t38 + t50 * t156 + t46 * t4 + ((-qJD(5) * t15 + t7) * t109 + t14 * t106 + t2 * t134) * t135 + ((qJD(5) * t14 + t6) * t109 + t15 * t106 + t176 * t134) * t138; 0, 0, 0, -t182, -t166, -qJ(2) * t166 + qJDD(2) - t223, -t131 * t166 - t173, t131 * t182 - t133 * t166, t33 * t131 + t32 * t133 + (-t110 * t132 + (t131 * t47 - t133 * t48) * t134) * qJD(1) - t222, 0, 0, 0, 0, 0, -t93 * t112 - t195 * t66 - t216 * t227, t94 * t112 - t195 * t69 - t215 * t227, 0, 0, 0, 0, 0, -(-t135 * t94 + t138 * t93) * t106 - t23 * t195 + (t135 * t161 + t138 * t162) * t109, (t135 * t93 + t138 * t94) * t106 - t156 * t195 + (-t135 * t162 + t138 * t161) * t109; 0, 0, 0, 0, 0, 0, (-t133 * t206 + t184) * t132, (qJDD(1) * t133 + t131 * t206) * t132, g(3) * t134 + ((t131 * t48 + t133 * t47) * qJD(1) - t160) * t132 + t90, 0, 0, 0, 0, 0, t35 + t210, t34 - t211, 0, 0, 0, 0, 0, t5 - t212, t4 + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t66, -t66 ^ 2 + t69 ^ 2, t34 + t211, -t35 + t210, -t112, -g(1) * t77 + g(2) * t75 + t121 * t218 - t155 * t227 - t85 * t69 + t147, g(1) * t78 - g(2) * t76 + t122 * t218 + t169 * t227 + t85 * t66 + t152, t229, t226, t230, t228, -t106, (-t10 * t135 - t213) * t109 + (-t106 * t138 + t109 * t190 - t23 * t69) * pkin(4) + t221, (-t10 * t109 - t176) * t138 + (t106 * t135 + t109 * t189 - t156 * t69) * pkin(4) + t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, t226, t230, t228, -t106, t109 * t158 + t221, (-t3 + (-qJD(5) - t109) * t8) * t138 + t225;];
tau_reg = t1;
