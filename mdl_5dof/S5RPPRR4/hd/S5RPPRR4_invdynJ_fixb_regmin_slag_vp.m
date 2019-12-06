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
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:44:59
% EndTime: 2019-12-05 17:45:10
% DurationCPUTime: 2.14s
% Computational Cost: add. (2038->276), mult. (5178->392), div. (0->0), fcn. (4071->14), ass. (0->172)
t135 = cos(pkin(8));
t190 = t135 * qJD(1);
t232 = qJD(4) - t190;
t110 = -qJD(5) - t232;
t136 = sin(qJ(5));
t139 = cos(qJ(5));
t133 = sin(pkin(8));
t132 = sin(pkin(9));
t134 = cos(pkin(9));
t137 = sin(qJ(4));
t140 = cos(qJ(4));
t94 = t140 * t132 + t137 * t134;
t153 = qJD(1) * t94;
t66 = t133 * t153;
t198 = qJD(1) * t133;
t181 = t132 * t198;
t168 = t137 * t181;
t205 = t140 * t134;
t184 = t133 * t205;
t69 = qJD(1) * t184 - t168;
t23 = t136 * t69 + t139 * t66;
t216 = t23 * t110;
t192 = qJD(5) * t139;
t193 = qJD(5) * t136;
t187 = qJDD(1) * t132;
t89 = t94 * qJD(4);
t34 = qJDD(1) * t184 + (-qJD(1) * t89 - t137 * t187) * t133;
t194 = qJD(4) * t140;
t180 = t134 * t194;
t35 = -qJD(4) * t168 + (qJD(1) * t180 + qJDD(1) * t94) * t133;
t4 = -t136 * t35 + t139 * t34 - t66 * t192 - t69 * t193;
t235 = t4 - t216;
t161 = -t136 * t66 + t139 * t69;
t234 = t161 * t23;
t217 = t161 * t110;
t5 = qJD(5) * t161 + t136 * t34 + t139 * t35;
t233 = -t5 - t217;
t231 = t161 ^ 2 - t23 ^ 2;
t150 = -t134 * t133 * pkin(6) + (-qJ(2) * t132 - pkin(3)) * t135;
t228 = -t135 * pkin(2) - pkin(1);
t97 = -t133 * qJ(3) + t228;
t84 = qJD(1) * t97 + qJD(2);
t74 = t134 * t84;
t36 = qJD(1) * t150 + t74;
t182 = qJ(2) * t190;
t48 = t132 * t84 + t134 * t182;
t41 = -pkin(6) * t181 + t48;
t160 = -t137 * t36 - t140 * t41;
t11 = -t66 * pkin(7) - t160;
t131 = pkin(9) + qJ(4);
t125 = qJ(5) + t131;
t120 = cos(t125);
t185 = t135 * qJDD(1);
t113 = -qJDD(4) + t185;
t188 = qJD(1) * qJD(2);
t178 = t135 * t188;
t196 = qJD(3) * t133;
t57 = -qJD(1) * t196 + qJDD(1) * t97 + qJDD(2);
t52 = t134 * t57;
t19 = qJDD(1) * t150 - t132 * t178 + t52;
t186 = t133 * qJDD(1);
t177 = t132 * t186;
t176 = t134 * t185;
t33 = qJ(2) * t176 + t132 * t57 + t134 * t178;
t21 = -pkin(6) * t177 + t33;
t174 = -t137 * t21 + t140 * t19;
t148 = qJD(4) * t160 + t174;
t2 = -t113 * pkin(4) - t34 * pkin(7) + t148;
t225 = g(1) * t133;
t111 = qJ(2) * t198 + qJD(3);
t85 = pkin(3) * t181 + t111;
t39 = t66 * pkin(4) + t85;
t119 = sin(t125);
t141 = cos(qJ(1));
t204 = t141 * t119;
t138 = sin(qJ(1));
t208 = t138 * t120;
t62 = t135 * t208 - t204;
t203 = t141 * t120;
t209 = t138 * t119;
t64 = -t135 * t203 - t209;
t9 = t11 * t193;
t230 = t39 * t23 + t120 * t225 + t9 + (t11 * t110 - t2) * t136 - g(2) * t62 - g(3) * t64;
t189 = qJ(2) * qJDD(1);
t218 = t139 * t11;
t173 = -t137 * t41 + t140 * t36;
t10 = -t69 * pkin(7) + t173;
t8 = pkin(4) * t232 + t10;
t162 = -t136 * t8 - t218;
t195 = qJD(4) * t137;
t155 = -t137 * t19 - t140 * t21 - t36 * t194 + t41 * t195;
t3 = -t35 * pkin(7) - t155;
t183 = -t136 * t3 + t139 * t2;
t61 = t135 * t209 + t203;
t63 = t135 * t204 - t208;
t227 = -g(2) * t61 + g(3) * t63 + qJD(5) * t162 + t119 * t225 - t39 * t161 + t183;
t224 = g(2) * t138;
t223 = g(3) * t141;
t92 = t134 * t97;
t44 = t92 + t150;
t212 = t132 * t133;
t60 = t134 * t135 * qJ(2) + t132 * t97;
t49 = -pkin(6) * t212 + t60;
t221 = t137 * t44 + t140 * t49;
t220 = -t135 * t153 + t89;
t93 = -t137 * t132 + t205;
t219 = t232 * t93;
t215 = t66 * t232;
t214 = t69 * t232;
t213 = qJDD(1) * pkin(1);
t211 = t132 * t135;
t142 = qJD(1) ^ 2;
t210 = t135 * t142;
t122 = sin(t131);
t207 = t138 * t122;
t123 = cos(t131);
t206 = t138 * t123;
t202 = t141 * t122;
t201 = t141 * t123;
t95 = pkin(3) * t212 + t133 * qJ(2);
t200 = -t132 ^ 2 - t134 ^ 2;
t128 = t133 ^ 2;
t199 = t135 ^ 2 + t128;
t197 = qJD(2) * t135;
t191 = t133 * qJD(2);
t179 = qJD(5) * t8 + t3;
t90 = qJ(2) * t186 + t133 * t188 + qJDD(3);
t172 = -t137 * t49 + t140 * t44;
t86 = -t132 * t197 - t134 * t196;
t87 = -t132 * t196 + t134 * t197;
t171 = -t137 * t87 + t140 * t86;
t170 = t199 * t142;
t169 = 0.2e1 * t199;
t58 = pkin(3) * t177 + t90;
t166 = qJD(5) * t94 + t220;
t165 = qJD(5) * t93 + t219;
t164 = g(2) * t141 + g(3) * t138;
t163 = -t223 + t224;
t82 = t94 * t133;
t83 = t93 * t133;
t37 = t136 * t83 + t139 * t82;
t38 = -t136 * t82 + t139 * t83;
t158 = t188 + t189;
t32 = -t158 * t211 + t52;
t59 = -qJ(2) * t211 + t92;
t157 = -t86 * qJD(1) - t59 * qJDD(1) - t32;
t156 = t87 * qJD(1) + t60 * qJDD(1) + t33;
t154 = t137 * t86 + t140 * t87 + t44 * t194 - t49 * t195;
t152 = -t164 - t213;
t121 = qJDD(2) - t213;
t151 = -t121 - t152;
t149 = t169 * t188 + t224;
t145 = t128 * t158 + t90 * t133 + t163;
t126 = t141 * qJ(2);
t107 = -qJDD(5) + t113;
t78 = -t135 * t201 - t207;
t77 = t135 * t202 - t206;
t76 = t135 * t206 - t202;
t75 = t135 * t207 + t201;
t72 = t133 * t180 - t195 * t212;
t71 = t133 * t89;
t50 = t72 * pkin(4) + t191;
t47 = -t132 * t182 + t74;
t46 = t82 * pkin(4) + t95;
t16 = t35 * pkin(4) + t58;
t15 = -t82 * pkin(7) + t221;
t14 = -t135 * pkin(4) - t83 * pkin(7) + t172;
t13 = qJD(5) * t38 - t136 * t71 + t139 * t72;
t12 = -qJD(5) * t37 - t136 * t72 - t139 * t71;
t7 = t71 * pkin(7) - qJD(4) * t221 + t171;
t6 = -t72 * pkin(7) + t154;
t1 = [qJDD(1), t164, -t163, t151 * t135, -t151 * t133, t169 * t189 + t149 - t223, -g(3) * t126 + (-t121 + t164) * pkin(1) + (t189 * t199 + t149) * qJ(2), (t134 * t164 + t157) * t135 + t145 * t132, (-t132 * t164 + t156) * t135 + t145 * t134, (-t132 * t156 + t134 * t157 + t164) * t133, t33 * t60 + t48 * t87 + t32 * t59 + t47 * t86 - g(2) * (-t138 * qJ(2) + t141 * t228) - g(3) * (t138 * t228 + t126) + (t90 * qJ(2) + qJ(3) * t164 + t111 * qJD(2)) * t133, t34 * t83 - t69 * t71, -t34 * t82 - t83 * t35 + t71 * t66 - t69 * t72, -t83 * t113 - t34 * t135 - t232 * t71, t82 * t113 + t35 * t135 - t232 * t72, t113 * t135, t171 * t232 - t172 * t113 - t174 * t135 + t66 * t191 + t95 * t35 + t58 * t82 + t85 * t72 - g(2) * t78 + g(3) * t76 + (-t135 * t160 - t221 * t232) * qJD(4), -g(2) * t77 - g(3) * t75 + t221 * t113 - t155 * t135 - t154 * t232 + t69 * t191 + t95 * t34 + t58 * t83 - t85 * t71, t12 * t161 + t4 * t38, -t12 * t23 - t13 * t161 - t4 * t37 - t38 * t5, -t38 * t107 - t12 * t110 - t4 * t135, t37 * t107 + t13 * t110 + t5 * t135, t107 * t135, -(-t136 * t6 + t139 * t7) * t110 - (-t136 * t15 + t139 * t14) * t107 - t183 * t135 + t50 * t23 + t46 * t5 + t16 * t37 + t39 * t13 - g(2) * t64 + g(3) * t62 + (-(-t136 * t14 - t139 * t15) * t110 - t162 * t135) * qJD(5), -g(2) * t63 - g(3) * t61 + t39 * t12 - t9 * t135 + t16 * t38 + t50 * t161 + t46 * t4 + ((-qJD(5) * t15 + t7) * t110 + t14 * t107 + t2 * t135) * t136 + ((qJD(5) * t14 + t6) * t110 + t15 * t107 + t179 * t135) * t139; 0, 0, 0, -t185, t186, -t170, -qJ(2) * t170 + qJDD(2) + t152, -t132 * t170 - t176, t132 * t185 - t134 * t170, t200 * t186, t33 * t132 + t32 * t134 + (-t111 * t133 + (t132 * t47 - t134 * t48) * t135) * qJD(1) - t164, 0, 0, 0, 0, 0, -t93 * t113 - t66 * t198 - t220 * t232, t94 * t113 - t69 * t198 - t219 * t232, 0, 0, 0, 0, 0, -(-t136 * t94 + t139 * t93) * t107 - t23 * t198 + (t136 * t165 + t139 * t166) * t110, (t136 * t93 + t139 * t94) * t107 - t161 * t198 + (-t136 * t166 + t139 * t165) * t110; 0, 0, 0, 0, 0, 0, 0, (-t134 * t210 + t187) * t133, (qJDD(1) * t134 + t132 * t210) * t133, t200 * t142 * t128, g(1) * t135 + ((t132 * t48 + t134 * t47) * qJD(1) + t163) * t133 + t90, 0, 0, 0, 0, 0, t35 + t214, t34 - t215, 0, 0, 0, 0, 0, t5 - t217, t4 + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t66, -t66 ^ 2 + t69 ^ 2, t34 + t215, -t35 + t214, -t113, -g(2) * t75 + g(3) * t77 + t122 * t225 - t160 * t232 - t85 * t69 + t148, -g(2) * t76 - g(3) * t78 + t123 * t225 + t173 * t232 + t85 * t66 + t155, t234, t231, t235, t233, -t107, (-t136 * t10 - t218) * t110 + (-t139 * t107 + t110 * t193 - t69 * t23) * pkin(4) + t227, (-t10 * t110 - t179) * t139 + (t136 * t107 + t110 * t192 - t161 * t69) * pkin(4) + t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, t231, t235, t233, -t107, t110 * t162 + t227, (-t3 + (-qJD(5) - t110) * t8) * t139 + t230;];
tau_reg = t1;
