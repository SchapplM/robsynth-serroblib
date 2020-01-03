% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:53
% EndTime: 2019-12-31 17:31:00
% DurationCPUTime: 2.61s
% Computational Cost: add. (2487->339), mult. (6446->504), div. (0->0), fcn. (5110->10), ass. (0->173)
t129 = cos(qJ(3));
t125 = sin(qJ(3));
t126 = sin(qJ(2));
t122 = sin(pkin(4));
t203 = qJD(1) * t122;
t181 = t126 * t203;
t163 = t125 * t181;
t123 = cos(pkin(4));
t202 = qJD(1) * t123;
t168 = qJD(2) + t202;
t64 = -t129 * t168 + t163;
t59 = qJD(4) + t64;
t190 = t126 * qJDD(1);
t174 = t122 * t190;
t130 = cos(qJ(2));
t193 = qJD(1) * qJD(2);
t175 = t130 * t193;
t240 = -t122 * t175 - t174;
t180 = t130 * t203;
t103 = -qJD(3) + t180;
t124 = sin(qJ(4));
t128 = cos(qJ(4));
t131 = cos(qJ(1));
t216 = t122 * t125;
t127 = sin(qJ(1));
t209 = t127 * t130;
t210 = t126 * t131;
t84 = t123 * t210 + t209;
t48 = t129 * t84 - t131 * t216;
t206 = t130 * t131;
t211 = t126 * t127;
t83 = -t123 * t206 + t211;
t239 = t124 * t48 - t128 * t83;
t238 = t124 * t83 + t128 * t48;
t192 = qJDD(1) * t123;
t114 = qJDD(2) + t192;
t148 = qJD(3) * t168;
t200 = qJD(2) * t130;
t178 = t125 * t200;
t196 = qJD(3) * t129;
t24 = -t129 * t114 + ((t126 * t196 + t178) * qJD(1) + t125 * t190) * t122 + t125 * t148;
t119 = t122 ^ 2;
t189 = 0.2e1 * t119;
t207 = t129 * t131;
t222 = t125 * t84;
t214 = t122 * t129;
t86 = -t123 * t211 + t206;
t50 = -t125 * t86 + t127 * t214;
t215 = t122 * t126;
t81 = -t123 * t129 + t125 * t215;
t144 = g(1) * t50 + g(2) * (-t122 * t207 - t222) - g(3) * t81;
t188 = pkin(1) * t202;
t78 = pkin(6) * t180 + t126 * t188;
t54 = t168 * pkin(7) + t78;
t149 = -pkin(2) * t130 - pkin(7) * t126 - pkin(1);
t73 = t149 * t122;
t58 = qJD(1) * t73;
t22 = t125 * t58 + t129 * t54;
t176 = t126 * t193;
t162 = t122 * t176;
t191 = qJDD(1) * t130;
t113 = t122 * t191;
t165 = qJD(2) * t188;
t187 = pkin(1) * t192;
t182 = -pkin(6) * t113 - t126 * t187 - t130 * t165;
t137 = -pkin(6) * t162 - t182;
t35 = pkin(7) * t114 + t137;
t159 = pkin(2) * t126 - pkin(7) * t130;
t145 = t159 * qJD(2);
t38 = (qJD(1) * t145 + t149 * qJDD(1)) * t122;
t135 = -t22 * qJD(3) - t125 * t35 + t129 * t38;
t74 = qJDD(3) - t113 + t162;
t4 = -pkin(3) * t74 - t135;
t66 = t125 * t168 + t129 * t181;
t237 = (pkin(3) * t66 + t59 * pkin(8)) * t59 + t144 + t4;
t20 = qJDD(4) + t24;
t99 = -pkin(3) * t129 - pkin(8) * t125 - pkin(2);
t236 = (t78 + t103 * (pkin(3) * t125 - pkin(8) * t129)) * t59 - t99 * t20;
t213 = t122 * t130;
t233 = pkin(1) * t126;
t205 = pkin(6) * t213 + t123 * t233;
t72 = pkin(7) * t123 + t205;
t153 = t125 * t73 + t129 * t72;
t77 = t122 * t145;
t115 = pkin(6) * t215;
t212 = t123 * t130;
t79 = (pkin(1) * t212 - t115) * qJD(2);
t234 = -t153 * qJD(3) - t125 * t79 + t129 * t77;
t150 = t103 * t124 - t128 * t66;
t23 = -qJD(3) * t163 + t125 * t114 + (t148 - t240) * t129;
t10 = -t150 * qJD(4) + t124 * t23 - t128 * t74;
t75 = -pkin(6) * t181 + t130 * t188;
t53 = -t168 * pkin(2) - t75;
t16 = t64 * pkin(3) - t66 * pkin(8) + t53;
t18 = -pkin(8) * t103 + t22;
t156 = t124 * t18 - t128 * t16;
t198 = qJD(3) * t125;
t147 = t125 * t38 + t129 * t35 + t58 * t196 - t54 * t198;
t3 = pkin(8) * t74 + t147;
t164 = t240 * pkin(6) - t126 * t165 + t130 * t187;
t36 = -pkin(2) * t114 - t164;
t6 = pkin(3) * t24 - pkin(8) * t23 + t36;
t1 = -t156 * qJD(4) + t124 * t6 + t128 * t3;
t132 = qJD(1) ^ 2;
t40 = t128 * t103 + t124 * t66;
t232 = t40 * t59;
t231 = t150 * t59;
t194 = qJD(4) * t128;
t195 = qJD(4) * t124;
t9 = -t103 * t194 + t124 * t74 + t128 * t23 - t66 * t195;
t230 = t9 * t124;
t76 = t159 * t203;
t227 = t125 * t76 + t129 * t75;
t226 = pkin(7) * qJD(3);
t225 = t103 * t64;
t224 = t124 * t20;
t221 = t128 * t20;
t219 = t66 * t103;
t218 = t103 * t125;
t217 = t119 * t132;
t208 = t129 * t130;
t120 = t126 ^ 2;
t204 = -t130 ^ 2 + t120;
t201 = qJD(2) * t126;
t199 = qJD(3) * t124;
t197 = qJD(3) * t128;
t186 = t59 * t199;
t185 = t59 * t197;
t184 = t130 * t217;
t183 = t124 * t213;
t179 = t122 * t201;
t177 = t122 * t123 * t132;
t169 = t128 * t59;
t167 = qJD(2) + 0.2e1 * t202;
t166 = t114 + t192;
t160 = -g(1) * t86 - g(2) * t84;
t8 = t124 * t16 + t128 * t18;
t71 = t115 + (-pkin(1) * t130 - pkin(2)) * t123;
t82 = t123 * t125 + t126 * t214;
t25 = pkin(3) * t81 - pkin(8) * t82 + t71;
t27 = -pkin(8) * t213 + t153;
t155 = -t124 * t27 + t128 * t25;
t154 = t124 * t25 + t128 * t27;
t21 = -t125 * t54 + t129 * t58;
t152 = -t125 * t72 + t129 * t73;
t45 = t124 * t82 + t128 * t213;
t146 = t125 * t77 + t129 * t79 + t73 * t196 - t72 * t198;
t85 = t123 * t209 + t210;
t141 = -g(1) * t85 - g(2) * t83 + g(3) * t213;
t80 = t205 * qJD(2);
t17 = pkin(3) * t103 - t21;
t139 = -pkin(8) * t20 + (t17 + t21) * t59;
t138 = -t141 - t36;
t136 = -pkin(7) * t74 - t103 * t53;
t2 = -t8 * qJD(4) - t124 * t3 + t128 * t6;
t134 = pkin(7) * qJD(4) * t59 + t141;
t133 = -g(3) * t215 + (pkin(8) * t181 - qJD(4) * t99 + t227) * t59 + t160;
t57 = (t124 * t126 + t128 * t208) * t203;
t56 = t124 * t129 * t180 - t128 * t181;
t51 = t127 * t216 + t129 * t86;
t46 = t128 * t82 - t183;
t44 = -t81 * qJD(3) + t200 * t214;
t43 = t82 * qJD(3) + t122 * t178;
t31 = t124 * t85 + t128 * t51;
t30 = -t124 * t51 + t128 * t85;
t28 = -pkin(3) * t181 + t125 * t75 - t129 * t76;
t26 = pkin(3) * t213 - t152;
t15 = -t45 * qJD(4) + t124 * t179 + t128 * t44;
t14 = -qJD(4) * t183 + t124 * t44 - t128 * t179 + t82 * t194;
t13 = pkin(3) * t43 - pkin(8) * t44 + t80;
t12 = -pkin(3) * t179 - t234;
t11 = pkin(8) * t179 + t146;
t5 = [qJDD(1), g(1) * t127 - g(2) * t131, g(1) * t131 + g(2) * t127, (qJDD(1) * t120 + 0.2e1 * t126 * t175) * t119, (t130 * t190 - t204 * t193) * t189, (t166 * t126 + t167 * t200) * t122, (t166 * t130 - t167 * t201) * t122, t114 * t123, -t80 * t168 - t115 * t114 + t164 * t123 + g(1) * t84 - g(2) * t86 + (t114 * t212 + (-t176 + t191) * t189) * pkin(1), -t79 * t168 - t205 * t114 - t137 * t123 - g(1) * t83 + g(2) * t85 + (-t175 - t190) * pkin(1) * t189, t23 * t82 + t44 * t66, -t23 * t81 - t24 * t82 - t43 * t66 - t44 * t64, -t103 * t44 + t74 * t82 + (-t130 * t23 + t66 * t201) * t122, t103 * t43 - t74 * t81 + (t130 * t24 - t201 * t64) * t122, (-t103 * t201 - t130 * t74) * t122, -t234 * t103 + t152 * t74 + t80 * t64 + t71 * t24 + t36 * t81 + t53 * t43 + g(1) * t48 - g(2) * t51 + (-t130 * t135 + t201 * t21) * t122, t146 * t103 - t153 * t74 + t80 * t66 + t71 * t23 + t36 * t82 + t53 * t44 - g(1) * t222 - g(2) * t50 + (-g(1) * t207 + t130 * t147 - t201 * t22) * t122, -t15 * t150 + t46 * t9, -t10 * t46 + t14 * t150 - t15 * t40 - t45 * t9, t15 * t59 - t150 * t43 + t20 * t46 + t81 * t9, -t10 * t81 - t14 * t59 - t20 * t45 - t40 * t43, t20 * t81 + t43 * t59, (-qJD(4) * t154 - t124 * t11 + t128 * t13) * t59 + t155 * t20 + t2 * t81 - t156 * t43 + t12 * t40 + t26 * t10 + t4 * t45 + t17 * t14 + g(1) * t238 - g(2) * t31, -(qJD(4) * t155 + t128 * t11 + t124 * t13) * t59 - t154 * t20 - t1 * t81 - t8 * t43 - t12 * t150 + t26 * t9 + t4 * t46 + t17 * t15 - g(1) * t239 - g(2) * t30; 0, 0, 0, -t126 * t184, t204 * t217, -t130 * t177 + t174, t126 * t177 + t113, t114, t78 * t168 + t217 * t233 - t141 + t164, t75 * t168 + pkin(1) * t184 + (pkin(6) * t193 + g(3)) * t215 - t160 + t182, t23 * t125 - t129 * t219, (t23 + t225) * t129 + (-t24 + t219) * t125, -t103 * t196 + t125 * t74 + (t103 * t208 - t126 * t66) * t203, t103 * t198 + t129 * t74 + (t126 * t64 - t130 * t218) * t203, t103 * t181, -t21 * t181 - pkin(2) * t24 - t78 * t64 + (-t75 * t103 + t136) * t125 + ((t76 + t226) * t103 + t138) * t129, -pkin(2) * t23 - t227 * t103 + t22 * t181 - t78 * t66 + t136 * t129 + (-t103 * t226 - t138) * t125, t125 * t128 * t9 - (-t125 * t195 + t128 * t196 - t57) * t150, t40 * t57 - t150 * t56 + (t124 * t150 - t128 * t40) * t196 + (-t10 * t128 - t230 + (t124 * t40 + t128 * t150) * qJD(4)) * t125, -t57 * t59 + (-t9 + t185) * t129 + (t103 * t150 - t195 * t59 + t221) * t125, t56 * t59 + (t10 - t186) * t129 + (t103 * t40 - t194 * t59 - t224) * t125, -t20 * t129 - t59 * t218, -t17 * t56 - t28 * t40 - t236 * t128 + t133 * t124 + (t17 * t199 - t2 + (qJD(3) * t40 - t224) * pkin(7) - t134 * t128) * t129 + (t17 * t194 + t4 * t124 + t103 * t156 + (t10 + t186) * pkin(7)) * t125, -t17 * t57 + t28 * t150 + t236 * t124 + t133 * t128 + (t17 * t197 + t1 + (-qJD(3) * t150 - t221) * pkin(7) + t134 * t124) * t129 + (-t17 * t195 + t4 * t128 + t103 * t8 + (t9 + t185) * pkin(7)) * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t64, -t64 ^ 2 + t66 ^ 2, t23 - t225, -t219 - t24, t74, -t103 * t22 - t53 * t66 + t135 - t144, g(1) * t51 + g(2) * t48 + g(3) * t82 - t103 * t21 + t53 * t64 - t147, -t150 * t169 + t230, (t9 - t232) * t128 + (-t10 + t231) * t124, t150 * t66 + t169 * t59 + t224, -t124 * t59 ^ 2 + t40 * t66 + t221, -t59 * t66, -pkin(3) * t10 + t139 * t124 - t237 * t128 + t156 * t66 - t22 * t40, -pkin(3) * t9 + t237 * t124 + t139 * t128 + t150 * t22 + t8 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150 * t40, t150 ^ 2 - t40 ^ 2, t9 + t232, -t10 - t231, t20, -g(1) * t30 + g(2) * t239 + g(3) * t45 + t17 * t150 + t8 * t59 + t2, g(1) * t31 + g(2) * t238 + g(3) * t46 - t156 * t59 + t17 * t40 - t1;];
tau_reg = t5;
