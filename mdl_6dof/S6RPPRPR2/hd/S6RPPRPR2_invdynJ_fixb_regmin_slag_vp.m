% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:40
% EndTime: 2019-03-09 01:42:45
% DurationCPUTime: 2.08s
% Computational Cost: add. (2397->315), mult. (5047->375), div. (0->0), fcn. (3742->14), ass. (0->175)
t124 = sin(pkin(10));
t126 = cos(pkin(10));
t130 = sin(qJ(4));
t223 = cos(qJ(4));
t93 = t223 * t124 + t130 * t126;
t237 = t93 * qJD(1);
t239 = qJD(6) + t237;
t129 = sin(qJ(6));
t241 = t129 * t239;
t132 = cos(qJ(6));
t197 = t130 * t124;
t180 = qJD(1) * t197;
t181 = t223 * t126;
t171 = qJD(1) * t181;
t176 = qJDD(1) * t223;
t184 = t126 * qJDD(1);
t183 = qJD(4) * t171 + t124 * t176 + t130 * t184;
t44 = qJD(4) * t180 - t183;
t43 = -qJDD(6) + t44;
t36 = t132 * t43;
t152 = -t239 * t241 - t36;
t189 = t129 * qJD(4);
t82 = -t171 + t180;
t63 = -t132 * t82 + t189;
t243 = t239 * t63;
t123 = qJ(1) + pkin(9);
t115 = sin(t123);
t117 = cos(t123);
t235 = -g(1) * t115 + g(2) * t117;
t127 = cos(pkin(9));
t218 = t127 * pkin(1);
t109 = -pkin(2) - t218;
t187 = qJDD(1) * t109;
t94 = qJDD(3) + t187;
t242 = -t235 - t94;
t122 = pkin(10) + qJ(4);
t114 = sin(t122);
t116 = cos(t122);
t170 = g(1) * t117 + g(2) * t115;
t142 = -g(3) * t114 - t170 * t116;
t179 = qJD(4) * t223;
t191 = qJD(4) * t130;
t111 = t126 * qJDD(2);
t125 = sin(pkin(9));
t104 = t125 * pkin(1) + qJ(3);
t90 = qJD(1) * qJD(3) + t104 * qJDD(1);
t55 = t111 + (-pkin(7) * qJDD(1) - t90) * t124;
t67 = t124 * qJDD(2) + t126 * t90;
t56 = pkin(7) * t184 + t67;
t113 = t126 * qJD(2);
t209 = pkin(7) * qJD(1);
t96 = t104 * qJD(1);
t61 = t113 + (-t96 - t209) * t124;
t75 = t124 * qJD(2) + t126 * t96;
t62 = t126 * t209 + t75;
t175 = -t130 * t55 - t61 * t179 + t62 * t191 - t223 * t56;
t240 = t142 - t175;
t210 = -t130 * t62 + t223 * t61;
t234 = qJD(5) - t210;
t185 = t124 * qJDD(1);
t160 = -t126 * t176 + t130 * t185;
t238 = 0.2e1 * t237 * qJD(4) + t160;
t236 = -qJD(6) + t239;
t24 = t130 * t61 + t223 * t62;
t22 = -qJD(4) * qJ(5) - t24;
t225 = pkin(5) * t82;
t12 = -t22 - t225;
t226 = pkin(4) + pkin(8);
t233 = t226 * t43 + (t12 - t24 + t225) * t239;
t105 = g(3) * t116;
t232 = -t170 * t114 + t105;
t212 = pkin(7) + t104;
t88 = t212 * t124;
t89 = t212 * t126;
t25 = (qJD(3) * t124 + qJD(4) * t89) * t130 - qJD(3) * t181 + t88 * t179;
t38 = -t130 * t88 + t223 * t89;
t231 = t25 * qJD(4) - t38 * qJDD(4) + t114 * t235;
t26 = t93 * qJD(3) + qJD(4) * t38;
t37 = t130 * t89 + t223 * t88;
t230 = -t26 * qJD(4) - t37 * qJDD(4) - t116 * t235;
t190 = qJD(6) * t132;
t92 = -t181 + t197;
t182 = t92 * t190;
t87 = t93 * qJD(4);
t229 = -(-t239 * t87 + t43 * t92) * t129 + t239 * t182;
t228 = t82 ^ 2;
t227 = t237 ^ 2;
t45 = qJD(1) * t87 + t160;
t224 = t45 * pkin(4);
t131 = sin(qJ(1));
t217 = t131 * pkin(1);
t108 = t126 * pkin(3) + pkin(2);
t95 = -t108 - t218;
t145 = -t93 * qJ(5) + t95;
t27 = t226 * t92 + t145;
t216 = t27 * t43;
t65 = t132 * qJD(4) + t129 * t82;
t215 = t65 * t82;
t214 = t82 * t63;
t213 = t237 * t82;
t17 = -qJD(6) * t189 + t132 * qJDD(4) + t129 * t45 + t82 * t190;
t86 = t124 * t191 - t126 * t179;
t211 = t17 * t93 - t65 * t86;
t207 = t132 * t239;
t206 = t17 * t132;
t205 = t82 * qJ(5);
t203 = qJD(6) * t92;
t202 = qJDD(4) * pkin(4);
t201 = t115 * t129;
t200 = t115 * t132;
t199 = t117 * t129;
t198 = t117 * t132;
t196 = t24 * qJD(4);
t194 = pkin(5) * t237 + t234;
t192 = t124 ^ 2 + t126 ^ 2;
t186 = qJDD(4) * qJ(5);
t11 = -t226 * qJD(4) + t194;
t78 = t95 * qJDD(1) + qJDD(3);
t140 = t44 * qJ(5) + t78;
t139 = -qJD(5) * t237 + t140;
t6 = t226 * t45 + t139;
t178 = qJD(6) * t11 + t6;
t81 = t95 * qJD(1) + qJD(3);
t141 = -qJ(5) * t237 + t81;
t19 = t226 * t82 + t141;
t150 = t130 * t56 + t62 * t179 + t61 * t191 - t223 * t55;
t146 = qJDD(5) + t150;
t2 = -t44 * pkin(5) - t226 * qJDD(4) + t146;
t177 = -qJD(6) * t19 + t2;
t173 = t129 * qJDD(4) - t132 * t45;
t133 = cos(qJ(1));
t168 = g(1) * t131 - g(2) * t133;
t18 = t65 * qJD(6) + t173;
t167 = -t18 * t93 + t63 * t86;
t165 = t237 * t87 - t44 * t92;
t164 = -t45 * t93 + t82 * t86;
t162 = pkin(4) * t116 + qJ(5) * t114;
t5 = t11 * t129 + t132 * t19;
t66 = -t124 * t90 + t111;
t159 = -t66 * t124 + t67 * t126;
t158 = t124 * (-t124 * t96 + t113) - t126 * t75;
t157 = t86 * qJ(5) - t93 * qJD(5);
t46 = -t86 * qJD(4) + t93 * qJDD(4);
t153 = t87 * qJD(4) + t92 * qJDD(4);
t151 = t108 + t162;
t149 = -t203 * t241 + t87 * t207 - t92 * t36;
t29 = t93 * pkin(5) + t37;
t7 = -qJD(4) * qJD(5) + t175 - t186;
t3 = -pkin(5) * t45 - t7;
t148 = t12 * t87 + t29 * t43 + t3 * t92;
t147 = t129 * t43 - t207 * t239;
t143 = -t187 + t242;
t138 = t3 + (t226 * t239 + t205) * t239 + t142;
t32 = t82 * pkin(4) + t141;
t137 = t237 * t32 + t146 + t232;
t128 = -pkin(7) - qJ(3);
t119 = t133 * pkin(1);
t77 = qJD(4) * t82;
t71 = -t114 * t201 + t198;
t70 = t114 * t200 + t199;
t69 = t114 * t199 + t200;
t68 = t114 * t198 - t201;
t42 = pkin(4) * t237 + t205;
t35 = t92 * pkin(4) + t145;
t33 = t87 * pkin(4) + t157;
t30 = -t92 * pkin(5) + t38;
t21 = -qJD(4) * pkin(4) + t234;
t20 = t226 * t87 + t157;
t16 = -t86 * pkin(5) + t26;
t15 = -t87 * pkin(5) - t25;
t9 = t139 + t224;
t8 = t146 - t202;
t4 = t11 * t132 - t129 * t19;
t1 = t132 * t2;
t10 = [qJDD(1), t168, g(1) * t133 + g(2) * t131 (t168 + (t125 ^ 2 + t127 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t143 * t126, -t143 * t124, t90 * t192 + t159 - t170, t94 * t109 - g(1) * (-pkin(2) * t115 + qJ(3) * t117 - t217) - g(2) * (pkin(2) * t117 + qJ(3) * t115 + t119) + t159 * t104 - t158 * qJD(3), -t237 * t86 - t44 * t93, t164 - t165, t46, -t153, 0, t95 * t45 + t78 * t92 + t81 * t87 + t230, -t95 * t44 + t78 * t93 - t81 * t86 + t231, -t21 * t86 + t22 * t87 + t237 * t26 + t25 * t82 - t37 * t44 - t38 * t45 + t7 * t92 + t8 * t93 - t170, -t32 * t87 - t33 * t82 - t35 * t45 - t9 * t92 - t230, -t237 * t33 + t32 * t86 + t35 * t44 - t9 * t93 - t231, g(1) * t217 - g(2) * t119 + t21 * t26 + t22 * t25 + t32 * t33 + t9 * t35 + t8 * t37 - t7 * t38 + (g(1) * t128 - g(2) * t151) * t117 + (g(1) * t151 + g(2) * t128) * t115, t65 * t182 + (t17 * t92 + t65 * t87) * t129 (-t129 * t63 + t132 * t65) * t87 + (-t129 * t18 + t206 + (-t129 * t65 - t132 * t63) * qJD(6)) * t92, t211 + t229, t149 + t167, -t239 * t86 - t43 * t93, -g(1) * t71 - g(2) * t69 + t1 * t93 + t15 * t63 + t30 * t18 - t4 * t86 + (-t20 * t239 - t6 * t93 + t216) * t129 + (t16 * t239 - t148) * t132 + ((-t129 * t29 - t132 * t27) * t239 - t5 * t93 + t12 * t129 * t92) * qJD(6), g(1) * t70 - g(2) * t68 + t15 * t65 + t30 * t17 + t5 * t86 + (-(qJD(6) * t29 + t20) * t239 + t216 - t178 * t93 + t12 * t203) * t132 + (-(-qJD(6) * t27 + t16) * t239 - t177 * t93 + t148) * t129; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, t124 * t67 + t126 * t66 - g(3), 0, 0, 0, 0, 0, -t153, -t46, t164 + t165, t153, t46, t21 * t87 + t22 * t86 - t7 * t93 + t8 * t92 - g(3), 0, 0, 0, 0, 0, t149 - t167, t211 - t229; 0, 0, 0, 0, -t184, t185, -t192 * qJD(1) ^ 2, t158 * qJD(1) - t242, 0, 0, 0, 0, 0, t238 (-t82 - t180) * qJD(4) + t183, -t227 - t228, -t238, t44 + t77, t224 - t22 * t82 + (-qJD(5) - t21) * t237 + t140 + t235, 0, 0, 0, 0, 0, t147 + t214, -t152 + t215; 0, 0, 0, 0, 0, 0, 0, 0, t213, t227 - t228 (t82 - t180) * qJD(4) + t183, -t160, qJDD(4), -t237 * t81 - t150 + t196 - t232, qJD(4) * t210 + t81 * t82 - t240, pkin(4) * t44 - qJ(5) * t45 + (-t22 - t24) * t237 + (t21 - t234) * t82, t42 * t82 + t137 - t196 - 0.2e1 * t202, 0.2e1 * t186 - t32 * t82 + t42 * t237 + (0.2e1 * qJD(5) - t210) * qJD(4) + t240, -t8 * pkin(4) - g(3) * t162 - t7 * qJ(5) - t234 * t22 - t21 * t24 - t32 * t42 + t170 * (pkin(4) * t114 - qJ(5) * t116) -t241 * t65 + t206 (-t239 * t65 - t18) * t132 + (-t17 + t243) * t129, t152 + t215, t147 - t214, t239 * t82, qJ(5) * t18 + t138 * t129 + t233 * t132 + t194 * t63 + t4 * t82, qJ(5) * t17 - t233 * t129 + t138 * t132 + t194 * t65 - t5 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 + t77, qJDD(4) - t213, -qJD(4) ^ 2 - t227, t22 * qJD(4) + t137 - t202, 0, 0, 0, 0, 0, -qJD(4) * t63 + t152, -qJD(4) * t65 + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t63, -t63 ^ 2 + t65 ^ 2, t17 + t243, t236 * t65 - t173, -t43, -g(1) * t68 - g(2) * t70 + t105 * t132 - t12 * t65 - t129 * t6 + t236 * t5 + t1, g(1) * t69 - g(2) * t71 + t12 * t63 + t4 * t239 - t178 * t132 + (-t177 - t105) * t129;];
tau_reg  = t10;
