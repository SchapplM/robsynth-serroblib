% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:32
% EndTime: 2019-12-05 16:13:43
% DurationCPUTime: 2.08s
% Computational Cost: add. (2613->214), mult. (5624->266), div. (0->0), fcn. (3609->8), ass. (0->155)
t128 = sin(qJ(3));
t130 = cos(qJ(3));
t167 = qJD(2) * t130;
t125 = sin(pkin(8));
t126 = cos(pkin(8));
t168 = qJD(2) * t128;
t94 = -t126 * qJD(3) + t125 * t168;
t162 = t94 * t167;
t166 = qJD(2) * qJD(3);
t157 = t130 * t166;
t165 = t128 * qJDD(2);
t101 = t157 + t165;
t75 = t125 * qJDD(3) + t126 * t101;
t137 = t75 + t162;
t113 = t128 * t166;
t164 = t130 * qJDD(2);
t102 = -t113 + t164;
t96 = t125 * qJD(3) + t126 * t168;
t194 = t96 * t94;
t142 = t102 - t194;
t182 = t126 * t142;
t121 = t130 ^ 2;
t133 = qJD(2) ^ 2;
t118 = t121 * t133;
t92 = t96 ^ 2;
t204 = -t92 - t118;
t40 = -t125 * t204 + t182;
t21 = t128 * t137 + t130 * t40;
t188 = t125 * t142;
t28 = -t126 * t204 - t188;
t242 = pkin(2) * t28 + pkin(6) * t21;
t129 = sin(qJ(2));
t131 = cos(qJ(2));
t241 = t129 * t21 + t131 * t28;
t239 = pkin(3) * t28;
t138 = -t75 + t162;
t153 = -t126 * qJDD(3) + t125 * t101;
t85 = t96 * t167;
t51 = t153 + t85;
t216 = -t125 * t138 - t126 * t51;
t196 = t94 ^ 2;
t49 = t92 + t196;
t228 = -t128 * t49 + t130 * t216;
t238 = pkin(6) * t228;
t237 = qJ(4) * t28;
t236 = qJ(4) * t40;
t218 = -t125 * t51 + t126 * t138;
t234 = t129 * t228 - t131 * t218;
t58 = t102 + t194;
t181 = t126 * t58;
t199 = -t196 - t118;
t206 = t125 * t199 - t181;
t187 = t125 * t58;
t205 = t126 * t199 + t187;
t52 = t153 - t85;
t214 = t128 * t52 + t130 * t205;
t232 = -pkin(2) * t206 + pkin(6) * t214;
t231 = pkin(3) * t49 + qJ(4) * t216;
t81 = t196 - t118;
t230 = t128 * (t126 * t81 + t188) + t130 * t51;
t229 = t129 * t214 - t131 * t206;
t225 = pkin(3) * t206;
t82 = -t92 + t118;
t223 = t126 * t82 - t187;
t222 = qJ(4) * t206;
t220 = -pkin(3) * t52 + qJ(4) * t205;
t219 = t125 * t81 - t182;
t190 = t125 * t137;
t215 = t128 * (t126 * t52 + t190) + t130 * (t92 - t196);
t213 = t128 * (-t125 * t82 - t181) + t130 * t138;
t212 = pkin(4) * t52;
t175 = sin(pkin(7));
t176 = cos(pkin(7));
t106 = -t176 * g(1) - t175 * g(2);
t122 = -g(3) + qJDD(1);
t78 = t131 * t106 + t129 * t122;
t67 = -t133 * pkin(2) + qJDD(2) * pkin(6) + t78;
t149 = -t130 * pkin(3) - t128 * qJ(4);
t99 = t149 * qJD(2);
t154 = qJD(2) * t99 + t67;
t201 = t154 * t128;
t148 = t101 + t157;
t77 = -t129 * t106 + t131 * t122;
t66 = -qJDD(2) * pkin(2) - t133 * pkin(6) - t77;
t36 = -t148 * qJ(4) + (-t102 + t113) * pkin(3) + t66;
t132 = qJD(3) ^ 2;
t140 = -t175 * g(1) + t176 * g(2);
t139 = t128 * t140;
t37 = -t132 * pkin(3) + qJDD(3) * qJ(4) + t154 * t130 + t139;
t193 = t125 * t36 + t126 * t37;
t64 = t94 * pkin(4) - t96 * qJ(5);
t200 = -t102 * qJ(5) - 0.2e1 * qJD(5) * t167 - t94 * t64 + t193;
t198 = -t125 * t52 + t126 * t137;
t98 = t130 * t140;
t151 = -qJDD(3) * pkin(3) - t132 * qJ(4) + qJDD(4) - t98;
t141 = t75 * qJ(5) - t151 - t212;
t177 = t130 * t94;
t179 = t128 * t67;
t197 = -(qJ(5) * t177 - t128 * t99) * qJD(2) - t141 + t179;
t195 = 2 * qJD(4);
t35 = t151 + t201;
t192 = t125 * t35;
t185 = t126 * t35;
t174 = qJ(5) * t125;
t173 = qJ(5) * t126;
t172 = qJD(4) * t94;
t171 = qJD(5) * t96;
t109 = t128 * t133 * t130;
t170 = t128 * (qJDD(3) + t109);
t169 = t130 * (qJDD(3) - t109);
t87 = -0.2e1 * t172;
t16 = t87 + t193;
t163 = t96 * t177;
t161 = pkin(4) * t126 + pkin(3);
t160 = t125 * t167;
t159 = t126 * t167;
t156 = t125 * t37 - t126 * t36;
t15 = t96 * t195 + t156;
t6 = t125 * t15 + t126 * t16;
t46 = -t98 + t179;
t47 = t130 * t67 + t139;
t23 = t128 * t46 + t130 * t47;
t152 = t94 * t160;
t150 = t128 * (t126 * t75 + t96 * t160) - t163;
t5 = t125 * t16 - t126 * t15;
t147 = t87 + t200;
t146 = -t126 * t153 - t152;
t76 = t96 * t159;
t145 = t76 + t152;
t144 = -pkin(2) + t149;
t89 = t130 * t102;
t135 = t89 + t128 * (-t125 * t96 + t126 * t94) * t167;
t10 = qJDD(5) - qJ(5) * t118 + t102 * pkin(4) + (t195 + t64) * t96 + t156;
t134 = t128 * (t125 * t153 - t94 * t159) + t163;
t120 = t128 ^ 2;
t117 = t120 * t133;
t105 = t117 + t118;
t104 = (t120 + t121) * qJDD(2);
t103 = -0.2e1 * t113 + t164;
t100 = 0.2e1 * t157 + t165;
t86 = 0.2e1 * t171;
t70 = -t169 - t128 * (-t117 - t132);
t69 = t130 * (-t118 - t132) - t170;
t43 = t125 * t75 - t76;
t13 = -0.2e1 * t171 + t197;
t12 = -t197 + t86 - t212;
t11 = t86 - t201 + (t137 + t162) * qJ(5) + t141;
t9 = -pkin(4) * t118 + t147;
t8 = qJ(5) * t49 + t10;
t7 = (t49 - t118) * pkin(4) + t147;
t4 = t128 * t35 + t130 * t6;
t3 = t125 * t10 + t126 * t9;
t2 = -t126 * t10 + t125 * t9;
t1 = t128 * t13 + t130 * t3;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t122, 0, 0, 0, 0, 0, 0, t131 * qJDD(2) - t129 * t133, -t129 * qJDD(2) - t131 * t133, 0, t129 * t78 + t131 * t77, 0, 0, 0, 0, 0, 0, t131 * t103 + t129 * t69, -t131 * t100 + t129 * t70, t129 * t104 + t131 * t105, t129 * t23 - t131 * t66, 0, 0, 0, 0, 0, 0, t229, t241, t234, t129 * t4 - t131 * t5, 0, 0, 0, 0, 0, 0, t229, t234, -t241, t129 * t1 - t131 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t77, -t78, 0, 0, t148 * t128, t130 * t100 + t128 * t103, t170 + t130 * (-t117 + t132), -t128 * t157 + t89, t128 * (t118 - t132) + t169, 0, pkin(2) * t103 + pkin(6) * t69 - t130 * t66, -pkin(2) * t100 + pkin(6) * t70 + t128 * t66, pkin(2) * t105 + pkin(6) * t104 + t23, -pkin(2) * t66 + pkin(6) * t23, t150, -t215, t213, t134, t230, t135, t128 * (t192 - t222) + t130 * (t15 - t225) + t232, t128 * (t185 + t237) + t130 * (t16 + t239) + t242, -t128 * t5 + t144 * t218 + t238, pkin(6) * t4 + t144 * t5, t150, t213, t215, t135, -t230, t134, t128 * (-t125 * t12 - t52 * t173 - t222) + t130 * (pkin(4) * t58 - qJ(5) * t199 + t10 - t225) + t232, t128 * (-qJ(4) * t218 - t125 * t7 + t126 * t8) + t130 * (-pkin(3) * t218 - pkin(4) * t138 + qJ(5) * t51) - pkin(2) * t218 + t238, t128 * (-pkin(4) * t190 + t126 * t11 - t237) + t130 * (-t239 + qJ(5) * t142 + 0.2e1 * t172 + (t204 + t118) * pkin(4) - t200) - t242, t128 * (-qJ(4) * t2 + (pkin(4) * t125 - t173) * t13) + t130 * (-pkin(3) * t2 + pkin(4) * t10 - qJ(5) * t9) - pkin(2) * t2 + pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, t117 - t118, t165, t109, t164, qJDD(3), -t46, -t47, 0, 0, t43, t198, t223, t146, t219, t145, -t185 + t220, -pkin(3) * t137 + t192 + t236, t231 + t6, -pkin(3) * t35 + qJ(4) * t6, t43, t223, -t198, t145, -t219, t146, t126 * t12 - t52 * t174 + t220, t125 * t8 + t126 * t7 + t231, t125 * t11 + t137 * t161 - t236, qJ(4) * t3 + (-t161 - t174) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t137, -t49, t35, 0, 0, 0, 0, 0, 0, t52, -t49, -t137, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t138, t204, t10;];
tauJ_reg = t14;
