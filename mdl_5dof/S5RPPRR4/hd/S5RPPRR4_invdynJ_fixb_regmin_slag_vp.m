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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:31:28
% EndTime: 2020-01-03 11:31:38
% DurationCPUTime: 2.25s
% Computational Cost: add. (2038->276), mult. (5178->395), div. (0->0), fcn. (4071->14), ass. (0->173)
t218 = qJDD(1) * pkin(1);
t121 = qJDD(2) - t218;
t142 = sin(qJ(1));
t145 = cos(qJ(1));
t204 = -g(2) * t145 - g(3) * t142;
t240 = t204 - t121;
t139 = cos(pkin(8));
t193 = t139 * qJD(1);
t236 = qJD(4) - t193;
t110 = -qJD(5) - t236;
t140 = sin(qJ(5));
t143 = cos(qJ(5));
t137 = sin(pkin(8));
t136 = sin(pkin(9));
t138 = cos(pkin(9));
t141 = sin(qJ(4));
t144 = cos(qJ(4));
t94 = t136 * t144 + t138 * t141;
t156 = qJD(1) * t94;
t66 = t137 * t156;
t201 = qJD(1) * t137;
t184 = t136 * t201;
t171 = t141 * t184;
t210 = t144 * t138;
t187 = t137 * t210;
t69 = qJD(1) * t187 - t171;
t23 = t140 * t69 + t143 * t66;
t221 = t23 * t110;
t195 = qJD(5) * t143;
t196 = qJD(5) * t140;
t190 = qJDD(1) * t136;
t89 = t94 * qJD(4);
t34 = qJDD(1) * t187 + (-qJD(1) * t89 - t141 * t190) * t137;
t197 = qJD(4) * t144;
t183 = t138 * t197;
t35 = -qJD(4) * t171 + (qJD(1) * t183 + qJDD(1) * t94) * t137;
t4 = -t140 * t35 + t143 * t34 - t195 * t66 - t196 * t69;
t239 = t4 - t221;
t164 = -t140 * t66 + t143 * t69;
t238 = t164 * t23;
t222 = t164 * t110;
t5 = qJD(5) * t164 + t140 * t34 + t143 * t35;
t237 = -t5 - t222;
t235 = t164 ^ 2 - t23 ^ 2;
t154 = -t138 * t137 * pkin(6) + (-qJ(2) * t136 - pkin(3)) * t139;
t227 = t139 * pkin(2);
t97 = -qJ(3) * t137 - pkin(1) - t227;
t84 = qJD(1) * t97 + qJD(2);
t74 = t138 * t84;
t36 = qJD(1) * t154 + t74;
t185 = qJ(2) * t193;
t48 = t136 * t84 + t138 * t185;
t41 = -pkin(6) * t184 + t48;
t163 = -t141 * t36 - t144 * t41;
t11 = -pkin(7) * t66 - t163;
t135 = pkin(9) + qJ(4);
t125 = qJ(5) + t135;
t120 = cos(t125);
t188 = t139 * qJDD(1);
t113 = -qJDD(4) + t188;
t191 = qJD(1) * qJD(2);
t181 = t139 * t191;
t199 = qJD(3) * t137;
t57 = -qJD(1) * t199 + qJDD(1) * t97 + qJDD(2);
t52 = t138 * t57;
t19 = qJDD(1) * t154 - t136 * t181 + t52;
t189 = t137 * qJDD(1);
t180 = t136 * t189;
t179 = t138 * t188;
t33 = qJ(2) * t179 + t136 * t57 + t138 * t181;
t21 = -pkin(6) * t180 + t33;
t177 = -t141 * t21 + t144 * t19;
t152 = qJD(4) * t163 + t177;
t2 = -t113 * pkin(4) - t34 * pkin(7) + t152;
t230 = g(1) * t137;
t111 = qJ(2) * t201 + qJD(3);
t85 = pkin(3) * t184 + t111;
t39 = pkin(4) * t66 + t85;
t119 = sin(t125);
t209 = t145 * t119;
t213 = t142 * t120;
t62 = t139 * t213 - t209;
t208 = t145 * t120;
t214 = t142 * t119;
t64 = t139 * t208 + t214;
t9 = t11 * t196;
t234 = t39 * t23 + t120 * t230 + t9 + g(2) * t62 + (t11 * t110 - t2) * t140 - g(3) * t64;
t192 = qJ(2) * qJDD(1);
t223 = t143 * t11;
t176 = -t141 * t41 + t144 * t36;
t10 = -pkin(7) * t69 + t176;
t8 = pkin(4) * t236 + t10;
t165 = -t140 * t8 - t223;
t198 = qJD(4) * t141;
t158 = -t141 * t19 - t144 * t21 - t197 * t36 + t198 * t41;
t3 = -pkin(7) * t35 - t158;
t186 = -t140 * t3 + t143 * t2;
t61 = -t139 * t214 - t208;
t63 = t139 * t209 - t213;
t232 = -g(2) * t61 - g(3) * t63 + qJD(5) * t165 + t119 * t230 - t164 * t39 + t186;
t229 = g(2) * t142;
t228 = g(3) * t145;
t92 = t138 * t97;
t44 = t92 + t154;
t217 = t136 * t137;
t60 = qJ(2) * t138 * t139 + t136 * t97;
t49 = -pkin(6) * t217 + t60;
t226 = t141 * t44 + t144 * t49;
t225 = -t139 * t156 + t89;
t93 = -t136 * t141 + t210;
t224 = t236 * t93;
t220 = t66 * t236;
t219 = t69 * t236;
t216 = t136 * t139;
t146 = qJD(1) ^ 2;
t215 = t139 * t146;
t122 = sin(t135);
t212 = t142 * t122;
t123 = cos(t135);
t211 = t142 * t123;
t207 = t145 * t122;
t206 = t145 * t123;
t95 = pkin(3) * t217 + qJ(2) * t137;
t205 = pkin(1) * t145 + qJ(2) * t142;
t203 = -t136 ^ 2 - t138 ^ 2;
t132 = t137 ^ 2;
t202 = t139 ^ 2 + t132;
t200 = qJD(2) * t139;
t194 = t137 * qJD(2);
t182 = qJD(5) * t8 + t3;
t90 = qJ(2) * t189 + t137 * t191 + qJDD(3);
t175 = -t141 * t49 + t144 * t44;
t86 = -t136 * t200 - t138 * t199;
t87 = -t136 * t199 + t138 * t200;
t174 = -t141 * t87 + t144 * t86;
t173 = t202 * t146;
t172 = 0.2e1 * t202;
t58 = pkin(3) * t180 + t90;
t169 = qJD(5) * t94 + t225;
t168 = qJD(5) * t93 + t224;
t166 = -t228 + t229;
t82 = t94 * t137;
t83 = t93 * t137;
t37 = t140 * t83 + t143 * t82;
t38 = -t140 * t82 + t143 * t83;
t161 = t191 + t192;
t32 = -t161 * t216 + t52;
t59 = -qJ(2) * t216 + t92;
t160 = -qJD(1) * t86 - qJDD(1) * t59 - t32;
t159 = t87 * qJD(1) + t60 * qJDD(1) + t33;
t157 = t141 * t86 + t144 * t87 + t197 * t44 - t198 * t49;
t155 = t218 + t240;
t153 = t172 * t191 + t228;
t149 = t132 * t161 + t90 * t137 - t166;
t127 = t142 * pkin(1);
t107 = -qJDD(5) + t113;
t78 = t139 * t206 + t212;
t77 = t139 * t207 - t211;
t76 = t139 * t211 - t207;
t75 = -t139 * t212 - t206;
t72 = t137 * t183 - t198 * t217;
t71 = t137 * t89;
t50 = pkin(4) * t72 + t194;
t47 = -t136 * t185 + t74;
t46 = pkin(4) * t82 + t95;
t16 = pkin(4) * t35 + t58;
t15 = -pkin(7) * t82 + t226;
t14 = -pkin(4) * t139 - pkin(7) * t83 + t175;
t13 = qJD(5) * t38 - t140 * t71 + t143 * t72;
t12 = -qJD(5) * t37 - t140 * t72 - t143 * t71;
t7 = t71 * pkin(7) - qJD(4) * t226 + t174;
t6 = -pkin(7) * t72 + t157;
t1 = [qJDD(1), t204, t166, t155 * t139, -t155 * t137, t172 * t192 + t153 - t229, -t121 * pkin(1) - g(2) * t205 - g(3) * t127 + (t192 * t202 + t153) * qJ(2), (t138 * t204 + t160) * t139 + t149 * t136, (-t136 * t204 + t159) * t139 + t149 * t138, (-t136 * t159 + t138 * t160 + t204) * t137, t33 * t60 + t48 * t87 + t32 * t59 + t47 * t86 - g(2) * (t145 * t227 + t205) - g(3) * (-t145 * qJ(2) + t142 * t227 + t127) + (t90 * qJ(2) + qJ(3) * t204 + t111 * qJD(2)) * t137, t34 * t83 - t69 * t71, -t34 * t82 - t35 * t83 + t66 * t71 - t69 * t72, -t113 * t83 - t139 * t34 - t236 * t71, t113 * t82 + t139 * t35 - t236 * t72, t113 * t139, t174 * t236 - t175 * t113 - t177 * t139 + t66 * t194 + t95 * t35 + t58 * t82 + t85 * t72 - g(2) * t78 - g(3) * t76 + (-t139 * t163 - t226 * t236) * qJD(4), g(2) * t77 - g(3) * t75 + t113 * t226 - t139 * t158 - t157 * t236 + t194 * t69 + t95 * t34 + t58 * t83 - t85 * t71, t12 * t164 + t38 * t4, -t12 * t23 - t13 * t164 - t37 * t4 - t38 * t5, -t107 * t38 - t110 * t12 - t139 * t4, t107 * t37 + t110 * t13 + t139 * t5, t107 * t139, -(-t140 * t6 + t143 * t7) * t110 - (t14 * t143 - t140 * t15) * t107 - t186 * t139 + t50 * t23 + t46 * t5 + t16 * t37 + t39 * t13 - g(2) * t64 - g(3) * t62 + (-(-t14 * t140 - t143 * t15) * t110 - t165 * t139) * qJD(5), g(2) * t63 - g(3) * t61 + t39 * t12 - t9 * t139 + t16 * t38 + t50 * t164 + t46 * t4 + ((-qJD(5) * t15 + t7) * t110 + t14 * t107 + t2 * t139) * t140 + ((qJD(5) * t14 + t6) * t110 + t15 * t107 + t182 * t139) * t143; 0, 0, 0, -t188, t189, -t173, -qJ(2) * t173 - t240, -t136 * t173 - t179, t136 * t188 - t138 * t173, t203 * t189, t33 * t136 + t32 * t138 + (-t111 * t137 + (t136 * t47 - t138 * t48) * t139) * qJD(1) - t204, 0, 0, 0, 0, 0, -t93 * t113 - t201 * t66 - t225 * t236, t94 * t113 - t201 * t69 - t224 * t236, 0, 0, 0, 0, 0, -(-t140 * t94 + t143 * t93) * t107 - t23 * t201 + (t140 * t168 + t143 * t169) * t110, (t140 * t93 + t143 * t94) * t107 - t164 * t201 + (-t140 * t169 + t143 * t168) * t110; 0, 0, 0, 0, 0, 0, 0, (-t138 * t215 + t190) * t137, (qJDD(1) * t138 + t136 * t215) * t137, t203 * t146 * t132, g(1) * t139 + ((t136 * t48 + t138 * t47) * qJD(1) - t166) * t137 + t90, 0, 0, 0, 0, 0, t35 + t219, t34 - t220, 0, 0, 0, 0, 0, t5 - t222, t4 + t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t66, -t66 ^ 2 + t69 ^ 2, t34 + t220, -t35 + t219, -t113, -g(2) * t75 - g(3) * t77 + t122 * t230 - t163 * t236 - t85 * t69 + t152, g(2) * t76 - g(3) * t78 + t123 * t230 + t176 * t236 + t85 * t66 + t158, t238, t235, t239, t237, -t107, (-t10 * t140 - t223) * t110 + (-t107 * t143 + t110 * t196 - t23 * t69) * pkin(4) + t232, (-t10 * t110 - t182) * t143 + (t107 * t140 + t110 * t195 - t164 * t69) * pkin(4) + t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, t235, t239, t237, -t107, t110 * t165 + t232, (-t3 + (-qJD(5) - t110) * t8) * t143 + t234;];
tau_reg = t1;
