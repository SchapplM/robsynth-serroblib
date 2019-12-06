% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:51
% EndTime: 2019-12-05 16:33:02
% DurationCPUTime: 2.75s
% Computational Cost: add. (2031->359), mult. (4912->537), div. (0->0), fcn. (3956->14), ass. (0->182)
t138 = cos(qJ(3));
t204 = qJD(2) * t138;
t118 = -qJD(5) + t204;
t134 = sin(qJ(5));
t137 = cos(qJ(5));
t130 = sin(pkin(10));
t132 = cos(pkin(10));
t198 = t132 * qJD(3);
t135 = sin(qJ(3));
t205 = qJD(2) * t135;
t94 = t130 * t205 - t198;
t203 = qJD(3) * t130;
t96 = t132 * t205 + t203;
t162 = t134 * t94 - t137 * t96;
t244 = t118 * t162;
t136 = sin(qJ(2));
t202 = qJD(3) * t135;
t192 = pkin(7) * t202;
t131 = sin(pkin(5));
t209 = qJD(1) * t131;
t139 = cos(qJ(2));
t212 = t138 * t139;
t168 = pkin(3) * t135 - qJ(4) * t138;
t76 = qJD(3) * t168 - qJD(4) * t135;
t232 = t130 * t192 + t132 * t76 - (-t130 * t212 + t132 * t136) * t209;
t243 = t130 * t76 - (t130 * t136 + t132 * t212) * t209;
t103 = qJD(2) * pkin(7) + t136 * t209;
t133 = cos(pkin(5));
t213 = t133 * t138;
t242 = qJD(1) * t213 - t135 * t103;
t126 = t138 * qJDD(2);
t196 = qJD(2) * qJD(3);
t180 = t135 * t196;
t241 = t180 - t126;
t240 = -t118 - qJD(5);
t239 = -qJDD(3) * pkin(3) + qJDD(4);
t140 = qJD(3) ^ 2;
t197 = qJD(1) * qJD(2);
t182 = t136 * t197;
t216 = t131 * t139;
t167 = -qJDD(1) * t216 + t131 * t182;
t224 = sin(pkin(9));
t173 = t224 * t136;
t225 = cos(pkin(9));
t174 = t225 * t139;
t78 = -t133 * t174 + t173;
t236 = g(2) * t78;
t172 = t224 * t139;
t175 = t225 * t136;
t80 = t133 * t172 + t175;
t237 = g(1) * t80;
t171 = t236 + t237;
t238 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t140 + (-g(3) * t139 + t182) * t131 - t167 + t171;
t122 = t132 * qJDD(3);
t181 = t138 * t196;
t193 = t135 * qJDD(2);
t149 = t181 + t193;
t63 = t130 * t149 - t122;
t194 = qJDD(3) * t130;
t64 = t132 * t149 + t194;
t9 = -qJD(5) * t162 + t134 * t64 + t137 * t63;
t235 = pkin(8) + qJ(4);
t195 = qJDD(1) * t133;
t179 = t135 * t195;
t70 = qJDD(2) * pkin(7) + (qJDD(1) * t136 + t139 * t197) * t131;
t14 = t179 + qJDD(3) * qJ(4) + t138 * t70 + (qJD(4) + t242) * qJD(3);
t159 = pkin(3) * t138 + qJ(4) * t135 + pkin(2);
t21 = qJD(2) * t76 - qJDD(2) * t159 + t167;
t7 = t130 * t21 + t132 * t14;
t214 = t132 * t138;
t160 = pkin(4) * t135 - pkin(8) * t214;
t234 = qJD(3) * t160 + t232;
t215 = t132 * t135;
t218 = t130 * t138;
t233 = (-pkin(7) * t215 - pkin(8) * t218) * qJD(3) + t243;
t208 = qJD(1) * t135;
t116 = t133 * t208;
t61 = t138 * t103 + t116;
t53 = qJD(3) * qJ(4) + t61;
t207 = qJD(1) * t139;
t189 = t131 * t207;
t62 = -qJD(2) * t159 - t189;
t18 = t130 * t62 + t132 * t53;
t231 = -t132 * t192 + t243;
t101 = t168 * qJD(2);
t27 = t130 * t101 + t132 * t242;
t219 = t130 * t134;
t99 = -t137 * t132 + t219;
t150 = t99 * t138;
t230 = qJD(2) * t150 - t99 * qJD(5);
t100 = t130 * t137 + t132 * t134;
t151 = t100 * t138;
t229 = -qJD(2) * t151 + t100 * qJD(5);
t228 = qJD(2) * pkin(2);
t226 = t134 * t96;
t32 = t137 * t94 + t226;
t227 = t118 * t32;
t66 = pkin(7) * t214 - t130 * t159;
t127 = pkin(10) + qJ(5);
t124 = sin(t127);
t221 = t124 * t138;
t125 = cos(t127);
t220 = t125 * t138;
t217 = t131 * t136;
t211 = qJDD(1) - g(3);
t128 = t135 ^ 2;
t210 = -t138 ^ 2 + t128;
t206 = qJD(2) * t131;
t201 = qJD(3) * t138;
t200 = qJD(5) * t135;
t199 = qJD(5) * t137;
t191 = pkin(4) * t130 + pkin(7);
t6 = -t130 * t14 + t132 * t21;
t2 = t241 * pkin(4) - pkin(8) * t64 + t6;
t5 = -pkin(8) * t63 + t7;
t190 = -t134 * t5 + t137 * t2;
t187 = t135 * t207;
t186 = t130 * t204;
t185 = t136 * t206;
t184 = t139 * t206;
t183 = qJ(4) * t126;
t17 = -t130 * t53 + t132 * t62;
t26 = t132 * t101 - t130 * t242;
t177 = t131 * t225;
t176 = t131 * t224;
t79 = t133 * t175 + t172;
t81 = -t133 * t173 + t174;
t170 = g(1) * t81 + g(2) * t79;
t169 = t134 * t2 + t137 * t5;
t11 = -pkin(4) * t204 - pkin(8) * t96 + t17;
t12 = -pkin(8) * t94 + t18;
t3 = t11 * t137 - t12 * t134;
t4 = t134 * t11 + t137 * t12;
t93 = t132 * t159;
t35 = -pkin(8) * t215 - t93 + (-pkin(7) * t130 - pkin(4)) * t138;
t48 = -pkin(8) * t130 * t135 + t66;
t166 = -t134 * t48 + t137 * t35;
t165 = t134 * t35 + t137 * t48;
t86 = t133 * t135 + t138 * t217;
t40 = -t130 * t86 - t132 * t216;
t41 = -t130 * t216 + t132 * t86;
t164 = -t134 * t41 + t137 * t40;
t163 = t134 * t40 + t137 * t41;
t141 = qJD(2) ^ 2;
t161 = qJDD(2) * t139 - t136 * t141;
t85 = t135 * t217 - t213;
t8 = -qJD(5) * t226 - t134 * t63 + t137 * t64 - t94 * t199;
t107 = t235 * t130;
t156 = pkin(8) * t186 + t132 * qJD(4) - qJD(5) * t107 - t27;
t108 = t235 * t132;
t155 = t160 * qJD(2) + t130 * qJD(4) + qJD(5) * t108 + t26;
t154 = g(1) * (t135 * t81 - t138 * t176) + g(2) * (t79 * t135 + t138 * t177) + g(3) * t85;
t43 = -t135 * t177 + t138 * t79;
t45 = t135 * t176 + t81 * t138;
t153 = g(1) * t45 + g(2) * t43 + g(3) * t86;
t152 = qJD(3) * t116 + t103 * t201 + t135 * t70 - t138 * t195;
t15 = t152 + t239;
t148 = -t15 + t154;
t147 = -g(3) * t217 - t170;
t146 = g(3) * t216 - t171;
t49 = -qJD(3) * pkin(3) + qJD(4) - t242;
t145 = -qJ(4) * t202 + (qJD(4) - t49) * t138;
t104 = -t189 - t228;
t143 = -pkin(7) * qJDD(3) + (t104 + t189 - t228) * qJD(3);
t142 = -t152 + t154;
t121 = -pkin(4) * t132 - pkin(3);
t102 = t191 * t135;
t98 = qJDD(5) + t241;
t89 = t191 * t201;
t74 = t99 * t135;
t73 = t100 * t135;
t65 = -pkin(7) * t218 - t93;
t47 = qJD(3) * t86 + t135 * t184;
t46 = -qJD(3) * t85 + t138 * t184;
t36 = pkin(4) * t186 + t61;
t30 = qJD(3) * t151 + t199 * t215 - t200 * t219;
t29 = -qJD(3) * t150 - t100 * t200;
t28 = pkin(4) * t94 + t49;
t25 = t130 * t185 + t132 * t46;
t24 = -t130 * t46 + t132 * t185;
t10 = pkin(4) * t63 + t15;
t1 = [t211, 0, t161 * t131, (-qJDD(2) * t136 - t139 * t141) * t131, 0, 0, 0, 0, 0, -qJD(3) * t47 - qJDD(3) * t85 + (t138 * t161 - t139 * t180) * t131, -qJD(3) * t46 - qJDD(3) * t86 + (-t135 * t161 - t139 * t181) * t131, -t40 * t126 + t47 * t94 + t63 * t85 + (-t138 * t24 + t40 * t202) * qJD(2), t41 * t126 + t47 * t96 + t64 * t85 + (t138 * t25 - t41 * t202) * qJD(2), -t24 * t96 - t25 * t94 - t40 * t64 - t41 * t63, t15 * t85 + t17 * t24 + t18 * t25 + t40 * t6 + t41 * t7 + t47 * t49 - g(3), 0, 0, 0, 0, 0, -(-qJD(5) * t163 - t134 * t25 + t137 * t24) * t118 + t164 * t98 + t47 * t32 + t85 * t9, -t47 * t162 + t85 * t8 + (qJD(5) * t164 + t134 * t24 + t137 * t25) * t118 - t163 * t98; 0, qJDD(2), t211 * t216 + t171, -t211 * t217 + t170, qJDD(2) * t128 + 0.2e1 * t138 * t180, 0.2e1 * t135 * t126 - 0.2e1 * t210 * t196, qJDD(3) * t135 + t138 * t140, qJDD(3) * t138 - t135 * t140, 0, t143 * t135 + t238 * t138, -t238 * t135 + t143 * t138, t147 * t130 + (-t94 * t189 + pkin(7) * t63 + t15 * t130 + (qJD(2) * t65 + t17) * qJD(3)) * t135 + (-t65 * qJDD(2) - t6 + (pkin(7) * t94 + t130 * t49) * qJD(3) - t232 * qJD(2) - t146 * t132) * t138, t147 * t132 + (-t96 * t189 + pkin(7) * t64 + t15 * t132 + (-qJD(2) * t66 - t18) * qJD(3)) * t135 + (t66 * qJDD(2) + t7 + (pkin(7) * t96 + t132 * t49) * qJD(3) + t231 * qJD(2) + t146 * t130) * t138, -t63 * t66 - t64 * t65 - t232 * t96 - t231 * t94 + (-t130 * t18 - t132 * t17) * t201 + (-t130 * t7 - t132 * t6 - t146) * t135, t6 * t65 + t7 * t66 + t231 * t18 + t232 * t17 + t159 * t237 + t159 * t236 + (t15 * t135 + t49 * t201 - t170) * pkin(7) + (-g(3) * pkin(7) * t136 + (-g(3) * t159 - t49 * t208) * t139) * t131, -t162 * t29 - t74 * t8, t162 * t30 - t29 * t32 - t73 * t8 + t74 * t9, -t118 * t29 - t138 * t8 - t162 * t202 - t74 * t98, t118 * t30 + t138 * t9 - t32 * t202 - t73 * t98, -t118 * t202 - t138 * t98, t166 * t98 - t190 * t138 + t3 * t202 + t89 * t32 + t102 * t9 + t10 * t73 + t28 * t30 - g(1) * (t124 * t81 - t80 * t220) - g(2) * (t124 * t79 - t78 * t220) + (t233 * t134 - t234 * t137) * t118 + (t118 * t165 + t138 * t4) * qJD(5) + (-t32 * t187 - g(3) * (t124 * t136 + t125 * t212)) * t131, -t89 * t162 + t102 * t8 - t10 * t74 + t28 * t29 - t165 * t98 + t169 * t138 - t4 * t202 - g(1) * (t125 * t81 + t80 * t221) - g(2) * (t125 * t79 + t78 * t221) + (t234 * t134 + t233 * t137) * t118 + (t118 * t166 + t138 * t3) * qJD(5) + (t162 * t187 - g(3) * (-t124 * t212 + t125 * t136)) * t131; 0, 0, 0, 0, -t135 * t141 * t138, t210 * t141, t193, t126, qJDD(3), qJD(3) * t61 - t104 * t205 + t142, -t179 + (-qJD(2) * t104 - t70) * t138 + t153, t130 * t183 - pkin(3) * t63 - t61 * t94 + t148 * t132 + (t130 * t145 - t17 * t135 + t138 * t26) * qJD(2), t132 * t183 - pkin(3) * t64 - t61 * t96 - t148 * t130 + (t132 * t145 + t18 * t135 - t138 * t27) * qJD(2), t26 * t96 + t27 * t94 + (-qJ(4) * t63 - qJD(4) * t94 + t17 * t204 + t7) * t132 + (qJ(4) * t64 + qJD(4) * t96 + t18 * t204 - t6) * t130 - t153, -t17 * t26 - t18 * t27 - t49 * t61 + (-t17 * t130 + t18 * t132) * qJD(4) + t148 * pkin(3) + (-t6 * t130 + t7 * t132 - t153) * qJ(4), t100 * t8 - t162 * t230, -t100 * t9 + t162 * t229 - t230 * t32 - t8 * t99, t100 * t98 - t230 * t118 + t162 * t205, t229 * t118 + t32 * t205 - t98 * t99, t118 * t205, (-t107 * t137 - t108 * t134) * t98 + t121 * t9 + t10 * t99 - t3 * t205 - t36 * t32 + t229 * t28 + (t134 * t156 + t137 * t155) * t118 + t154 * t125, t121 * t8 + t10 * t100 - (-t107 * t134 + t108 * t137) * t98 + t36 * t162 + t4 * t205 + t230 * t28 + (-t134 * t155 + t137 * t156) * t118 - t154 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130 * t193 - t122 + (-t96 + t203) * t204, t132 * t193 + t194 + (t94 + t198) * t204, -t94 ^ 2 - t96 ^ 2, t17 * t96 + t18 * t94 - t142 + t239, 0, 0, 0, 0, 0, t9 + t244, t8 + t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162 * t32, t162 ^ 2 - t32 ^ 2, t8 - t227, -t9 + t244, t98, t28 * t162 - g(1) * (-t124 * t45 + t125 * t80) - g(2) * (-t124 * t43 + t125 * t78) - g(3) * (-t124 * t86 - t125 * t216) + t190 + t240 * t4, t28 * t32 - g(1) * (-t124 * t80 - t125 * t45) - g(2) * (-t124 * t78 - t125 * t43) - g(3) * (t124 * t216 - t125 * t86) - t169 + t240 * t3;];
tau_reg = t1;
