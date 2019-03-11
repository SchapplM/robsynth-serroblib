% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:13:57
% EndTime: 2019-03-09 02:14:03
% DurationCPUTime: 2.39s
% Computational Cost: add. (3501->341), mult. (7019->425), div. (0->0), fcn. (4961->10), ass. (0->174)
t120 = pkin(9) + qJ(4);
t107 = sin(t120);
t108 = cos(t120);
t131 = sin(qJ(1));
t134 = cos(qJ(1));
t231 = g(1) * t131 - g(2) * t134;
t139 = -g(3) * t107 + t231 * t108;
t130 = sin(qJ(4));
t133 = cos(qJ(4));
t188 = qJD(4) * t133;
t189 = qJD(4) * t130;
t124 = sin(pkin(9));
t128 = -pkin(1) - qJ(3);
t226 = -qJD(1) * qJD(3) + qJDD(1) * t128;
t87 = qJDD(2) + t226;
t170 = -pkin(7) * qJDD(1) + t87;
t63 = t170 * t124;
t125 = cos(pkin(9));
t64 = t170 * t125;
t95 = t128 * qJD(1) + qJD(2);
t173 = -pkin(7) * qJD(1) + t95;
t68 = t173 * t124;
t69 = t173 * t125;
t152 = -t130 * t63 + t133 * t64 - t68 * t188 - t69 * t189;
t14 = -qJDD(4) * pkin(4) - t152;
t86 = t124 * t133 + t125 * t130;
t232 = t86 * qJD(1);
t240 = qJD(5) + t232;
t243 = pkin(8) * qJD(5) * t240 + t139 + t14;
t129 = sin(qJ(5));
t132 = cos(qJ(5));
t196 = t133 * t125;
t85 = t124 * t130 - t196;
t143 = t85 * qJDD(1);
t186 = qJD(5) * t132;
t190 = qJD(4) * t129;
t191 = qJD(1) * t124;
t176 = t130 * t191;
t178 = qJD(1) * t196;
t80 = -t176 + t178;
t24 = -t132 * qJDD(4) - t129 * t143 + t80 * t186 + (qJD(5) - t232) * t190;
t167 = t132 * t240;
t141 = -qJD(4) * t176 + t86 * qJDD(1);
t177 = t125 * t188;
t50 = qJD(1) * t177 + t141;
t47 = qJDD(5) + t50;
t206 = t129 * t47;
t242 = -t167 * t240 - t206;
t192 = t124 ^ 2 + t125 ^ 2;
t238 = t192 * t95;
t236 = -t130 * t68 + t133 * t69;
t162 = g(1) * t134 + g(2) * t131;
t121 = qJDD(1) * qJ(2);
t122 = qJD(1) * qJD(2);
t230 = t121 + t122;
t93 = qJDD(3) + t230;
t235 = t93 - t162;
t234 = qJD(4) * t232;
t233 = t162 * t108;
t197 = t132 * t134;
t200 = t129 * t131;
t72 = -t107 * t200 + t197;
t198 = t131 * t132;
t199 = t129 * t134;
t74 = t107 * t199 + t198;
t229 = -g(1) * t72 - g(2) * t74;
t155 = t130 * t64 + t133 * t63;
t13 = qJDD(4) * pkin(8) + qJD(4) * t236 + t155;
t35 = t130 * t69 + t133 * t68;
t32 = qJD(4) * pkin(8) + t35;
t106 = qJD(1) * qJ(2) + qJD(3);
t90 = pkin(3) * t191 + t106;
t33 = pkin(4) * t232 - pkin(8) * t80 + t90;
t16 = t129 * t33 + t132 * t32;
t136 = -t143 - t234;
t183 = t124 * qJDD(1);
t84 = pkin(3) * t183 + t93;
t20 = t50 * pkin(4) - t136 * pkin(8) + t84;
t19 = t132 * t20;
t185 = t132 * qJD(4);
t187 = qJD(5) * t129;
t23 = -qJD(5) * t185 - t129 * qJDD(4) - t132 * t136 + t80 * t187;
t59 = t132 * t80 + t190;
t1 = pkin(5) * t47 + qJ(6) * t23 - t16 * qJD(5) - qJD(6) * t59 - t129 * t13 + t19;
t57 = t129 * t80 - t185;
t8 = -qJ(6) * t57 + t16;
t228 = t240 * t8 + t1;
t217 = g(3) * t108;
t227 = t107 * t231 + t217;
t225 = t59 ^ 2;
t224 = 0.2e1 * t122;
t15 = -t129 * t32 + t132 * t33;
t7 = -qJ(6) * t59 + t15;
t6 = pkin(5) * t240 + t7;
t223 = -t7 + t6;
t126 = -qJ(6) - pkin(8);
t172 = qJD(5) * t126;
t201 = qJ(6) * t132;
t49 = pkin(4) * t80 + pkin(8) * t232;
t43 = t132 * t49;
t220 = -pkin(5) * t80 + t132 * t172 - t232 * t201 - t43 + (-qJD(6) + t236) * t129;
t219 = pkin(5) * t129;
t112 = t124 * pkin(3);
t216 = t57 * t232;
t215 = t57 * t80;
t214 = t59 * t80;
t213 = -pkin(7) + t128;
t212 = -t129 * t24 - t57 * t186;
t211 = t129 * t49 + t132 * t236;
t100 = qJ(2) + t112;
t48 = pkin(4) * t86 + pkin(8) * t85 + t100;
t88 = t213 * t124;
t89 = t213 * t125;
t53 = t130 * t89 + t133 * t88;
t51 = t132 * t53;
t210 = t129 * t48 + t51;
t81 = -t124 * t188 - t125 * t189;
t209 = t81 * qJD(4) - t85 * qJDD(4);
t202 = qJ(6) * t129;
t208 = qJD(6) * t132 + t129 * t172 - t202 * t232 - t211;
t207 = t129 * t23;
t205 = t129 * t59;
t40 = t132 * t47;
t203 = pkin(1) * qJDD(1);
t194 = t134 * pkin(1) + t131 * qJ(2);
t52 = t130 * t88 - t133 * t89;
t28 = -t86 * qJD(3) - t52 * qJD(4);
t82 = -t124 * t189 + t177;
t45 = pkin(4) * t82 - pkin(8) * t81 + qJD(2);
t182 = t129 * t45 + t132 * t28 + t48 * t186;
t180 = t85 * t187;
t179 = t85 * t186;
t175 = g(2) * t194;
t174 = pkin(7) + qJ(3) + t219;
t171 = t192 * t87;
t169 = -qJD(5) * t33 - t13;
t168 = t129 * t240;
t166 = qJD(5) * t86 + qJD(1);
t165 = qJDD(2) - t203;
t164 = -t32 * t186 + t19;
t146 = t129 * t20 + t132 * t13 + t33 * t186 - t32 * t187;
t2 = -qJ(6) * t24 - qJD(6) * t57 + t146;
t163 = -t240 * t6 + t2;
t31 = -qJD(4) * pkin(4) - t236;
t160 = -t14 * t85 + t31 * t81;
t159 = t23 * t85 + t59 * t81;
t158 = -t23 * t86 + t59 * t82;
t157 = -t24 * t86 - t57 * t82;
t156 = -t240 * t81 + t47 * t85;
t153 = -qJ(6) * t81 + qJD(6) * t85;
t105 = pkin(5) * t132 + pkin(4);
t150 = -t105 * t107 - t108 * t126;
t149 = -qJD(4) * t82 - qJDD(4) * t86;
t148 = t40 + (-t129 * t232 - t187) * t240;
t147 = -pkin(8) * t47 + t240 * t31;
t142 = -t150 + t112;
t5 = pkin(5) * t24 + qJDD(6) + t14;
t138 = t230 + t235;
t29 = -t85 * qJD(3) + t53 * qJD(4);
t135 = qJD(1) ^ 2;
t114 = t134 * qJ(2);
t92 = t126 * t132;
t91 = t126 * t129;
t75 = t107 * t197 - t200;
t73 = t107 * t198 + t199;
t56 = t57 ^ 2;
t41 = t132 * t48;
t38 = t132 * t45;
t21 = pkin(5) * t57 + qJD(6) + t31;
t17 = t85 * t202 + t210;
t11 = pkin(5) * t86 - t129 * t53 + t85 * t201 + t41;
t4 = qJ(6) * t179 + (-qJD(5) * t53 + t153) * t129 + t182;
t3 = pkin(5) * t82 - t129 * t28 + t38 + t153 * t132 + (-t51 + (-qJ(6) * t85 - t48) * t129) * qJD(5);
t9 = [qJDD(1), t231, t162, qJDD(2) - t231 - 0.2e1 * t203, 0.2e1 * t121 + t224 - t162, -t165 * pkin(1) - g(1) * (-pkin(1) * t131 + t114) - t175 + (t121 + t224) * qJ(2), t138 * t124, t138 * t125, t231 + t192 * (-t226 - t87) t93 * qJ(2) + t106 * qJD(2) - g(1) * (t128 * t131 + t114) - g(2) * (qJ(3) * t134 + t194) + t128 * t171 - qJD(3) * t238, -t136 * t85 + t80 * t81, -t136 * t86 - t232 * t81 + t85 * t50 - t80 * t82, t209, t149, 0, qJD(2) * t232 - qJD(4) * t29 - qJDD(4) * t52 + t100 * t50 - t162 * t107 + t82 * t90 + t84 * t86, qJD(2) * t80 - t53 * qJDD(4) + t90 * t81 - t84 * t85 - t233 - t100 * t143 + (-t100 * t232 - t28) * qJD(4), t159 * t132 + t59 * t180 (-t132 * t57 - t205) * t81 + (-t207 + t132 * t24 + (-t129 * t57 + t132 * t59) * qJD(5)) * t85, -t156 * t132 + t180 * t240 + t158, t129 * t156 + t179 * t240 + t157, t240 * t82 + t47 * t86 (-t186 * t53 + t38) * t240 + t41 * t47 + t164 * t86 + t15 * t82 + t29 * t57 + t52 * t24 - t31 * t179 - g(1) * t75 - g(2) * t73 + ((-qJD(5) * t48 - t28) * t240 - t53 * t47 + t169 * t86 + t160) * t129 -(-t187 * t53 + t182) * t240 - t210 * t47 - t146 * t86 - t16 * t82 + t29 * t59 - t52 * t23 + t31 * t180 + g(1) * t74 - g(2) * t72 + t160 * t132, t11 * t23 - t17 * t24 - t3 * t59 - t4 * t57 + (-t129 * t8 - t132 * t6) * t81 + t233 + (t1 * t132 + t129 * t2 + (-t129 * t6 + t132 * t8) * qJD(5)) * t85, t2 * t17 + t8 * t4 + t1 * t11 + t6 * t3 + t5 * (-t219 * t85 + t52) + t21 * ((t129 * t81 - t179) * pkin(5) + t29) - g(1) * t114 - t175 + (-g(1) * t142 - g(2) * t174) * t134 + (-g(1) * (-pkin(1) - t174) - g(2) * t142) * t131; 0, 0, 0, qJDD(1), -t135, -qJ(2) * t135 + t165 - t231, -t135 * t124, -t135 * t125, -t192 * qJDD(1), -qJD(1) * t106 + t171 - t231, 0, 0, 0, 0, 0, -qJD(1) * t232 + t209, -qJD(1) * t80 + t149, 0, 0, 0, 0, 0, -t86 * t206 + t24 * t85 - t57 * t81 + (-t129 * t82 - t132 * t166) * t240, -t86 * t40 + (t129 * t166 - t132 * t82) * t240 - t159 (t166 * t59 + t157) * t132 + (t166 * t57 + t158) * t129, -t21 * t81 + t5 * t85 + (-t166 * t6 + t2 * t86 + t8 * t82) * t132 + (-t1 * t86 - t166 * t8 - t6 * t82) * t129 - t231; 0, 0, 0, 0, 0, 0, t183, t125 * qJDD(1), -t192 * t135, qJD(1) * t238 + t235, 0, 0, 0, 0, 0 (t80 + t178) * qJD(4) + t141, -t143 - 0.2e1 * t234, 0, 0, 0, 0, 0, t148 - t215, -t214 + t242 (t23 - t216) * t132 + t59 * t168 + t212, t163 * t129 + t228 * t132 - t21 * t80 - t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t232, -t232 ^ 2 + t80 ^ 2, -t143 (t80 - t178) * qJD(4) - t141, qJDD(4), qJD(4) * t35 - t80 * t90 - t139 + t152, t232 * t90 - t155 + t227, t59 * t167 - t207 (-t23 - t216) * t132 - t240 * t205 + t212, -t214 - t242, t148 + t215, -t240 * t80, -pkin(4) * t24 - t15 * t80 - t35 * t57 - t43 * t240 + (t236 * t240 + t147) * t129 - t243 * t132, pkin(4) * t23 + t243 * t129 + t147 * t132 + t16 * t80 + t211 * t240 - t35 * t59, -t228 * t129 + t163 * t132 - t208 * t57 - t220 * t59 + t23 * t91 + t24 * t92 - t227, -t2 * t92 + t1 * t91 - t5 * t105 - g(3) * t150 + t208 * t8 + t220 * t6 + (pkin(5) * t168 - t35) * t21 - t231 * (t105 * t108 - t107 * t126); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t57, -t56 + t225, t240 * t57 - t23, t240 * t59 - t24, t47, t16 * t240 - t31 * t59 + (t169 + t217) * t129 + t164 + t229, g(1) * t73 - g(2) * t75 + t132 * t217 + t15 * t240 + t31 * t57 - t146, pkin(5) * t23 - t223 * t57, t223 * t8 + (t129 * t217 - t21 * t59 + t1 + t229) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56 - t225, t57 * t8 + t59 * t6 + t139 + t5;];
tau_reg  = t9;
