% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRR4
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:41
% EndTime: 2019-12-05 15:51:46
% DurationCPUTime: 2.57s
% Computational Cost: add. (2677->338), mult. (6484->488), div. (0->0), fcn. (5398->12), ass. (0->173)
t123 = sin(qJ(4));
t126 = cos(qJ(4));
t160 = t126 * pkin(4) + t123 * pkin(8);
t151 = -pkin(3) - t160;
t116 = sin(pkin(10));
t119 = cos(pkin(10));
t124 = sin(qJ(2));
t127 = cos(qJ(2));
t118 = sin(pkin(5));
t196 = qJD(2) * t118;
t172 = qJD(1) * t196;
t185 = qJDD(1) * t118;
t230 = t124 * t185 + t127 * t172;
t213 = qJDD(2) * pkin(2);
t96 = t127 * t185;
t64 = -t124 * t172 + t213 + t96;
t30 = -t230 * t116 + t119 * t64;
t159 = pkin(4) * t123 - pkin(8) * t126;
t82 = t159 * qJD(4);
t12 = qJD(2) * t82 + qJDD(2) * t151 - t30;
t122 = sin(qJ(5));
t125 = cos(qJ(5));
t121 = cos(pkin(5));
t100 = t121 * qJD(1) + qJD(3);
t197 = qJD(1) * t118;
t178 = t124 * t197;
t174 = t127 * t197;
t84 = qJD(2) * pkin(2) + t174;
t55 = t116 * t84 + t119 * t178;
t53 = qJD(2) * pkin(7) + t55;
t35 = t123 * t100 + t126 * t53;
t33 = qJD(4) * pkin(8) + t35;
t89 = t116 * t178;
t54 = t119 * t84 - t89;
t39 = qJD(2) * t151 - t54;
t155 = t122 * t33 - t125 * t39;
t193 = qJD(4) * t126;
t31 = t116 * t64 + t230 * t119;
t29 = qJDD(2) * pkin(7) + t31;
t97 = t121 * qJDD(1) + qJDD(3);
t183 = -t100 * t193 - t123 * t97 - t126 * t29;
t194 = qJD(4) * t123;
t7 = -t194 * t53 - t183;
t5 = qJDD(4) * pkin(8) + t7;
t1 = -t155 * qJD(5) + t122 * t12 + t125 * t5;
t187 = t126 * qJD(2);
t102 = -qJD(5) + t187;
t232 = -t155 * t102 + t1;
t152 = t127 * t116 + t124 * t119;
t68 = t152 * t118;
t60 = qJD(1) * t68;
t231 = -t60 + t82;
t191 = qJD(5) * t123;
t229 = qJD(2) * t191 - qJDD(4);
t184 = t123 * qJDD(2);
t47 = (qJD(4) * (qJD(5) + t187) + t184) * t122 + t229 * t125;
t10 = t122 * t39 + t125 * t33;
t2 = -qJD(5) * t10 + t125 * t12 - t122 * t5;
t228 = -t10 * t102 + t2;
t74 = t124 * t116 - t127 * t119;
t88 = t126 * t97;
t8 = -t35 * qJD(4) - t123 * t29 + t88;
t188 = t125 * qJD(4);
t173 = t126 * t188;
t205 = t123 * t125;
t111 = t126 * qJDD(2);
t186 = qJD(2) * qJD(4);
t169 = t123 * t186;
t73 = qJDD(5) - t111 + t169;
t227 = -t102 * (-t122 * t191 + t173) + t73 * t205;
t225 = t119 * pkin(2);
t195 = qJD(2) * t123;
t76 = t122 * t195 - t188;
t189 = t122 * qJD(4);
t78 = t125 * t195 + t189;
t224 = t78 * t76;
t105 = t116 * pkin(2) + pkin(7);
t176 = t105 * t194;
t203 = t125 * t126;
t206 = t122 * t126;
t72 = t151 - t225;
t48 = -t105 * t206 + t125 * t72;
t63 = t119 * t174 - t89;
t222 = qJD(5) * t48 + t231 * t122 - t125 * t176 - t63 * t203;
t49 = t105 * t203 + t122 * t72;
t221 = -qJD(5) * t49 + t122 * t176 + t231 * t125 + t63 * t206;
t220 = -t76 * t173 - t47 * t205;
t218 = t123 * t53;
t217 = t76 * t102;
t216 = t78 * t102;
t215 = qJD(4) * t76;
t214 = qJD(5) * t76;
t117 = sin(pkin(9));
t212 = t117 * t124;
t211 = t118 * t123;
t210 = t118 * t126;
t209 = t118 * t127;
t208 = t121 * t124;
t207 = t121 * t127;
t201 = t63 * qJD(2);
t200 = qJDD(1) - g(3);
t114 = t123 ^ 2;
t115 = t126 ^ 2;
t199 = t114 - t115;
t192 = qJD(5) * t122;
t190 = qJD(5) * t125;
t181 = t78 * t193;
t120 = cos(pkin(9));
t180 = t120 * t207;
t129 = qJD(2) ^ 2;
t179 = t123 * t129 * t126;
t177 = t102 * t189;
t175 = t123 * t190;
t46 = -qJD(2) * t173 - qJD(5) * t188 + t229 * t122 - t125 * t184;
t167 = t46 * t126 + t78 * t194;
t166 = -t46 + t214;
t164 = t78 * t175;
t67 = t74 * t118;
t163 = pkin(2) * t209 - t67 * pkin(3) + t68 * pkin(7);
t161 = t126 * t169;
t158 = t10 * t125 + t122 * t155;
t157 = -t10 * t122 + t125 * t155;
t69 = t152 * t121;
t40 = t117 * t74 - t120 * t69;
t43 = t117 * t69 + t120 * t74;
t51 = t121 * t123 + t68 * t126;
t20 = t67 * t122 + t51 * t125;
t19 = -t51 * t122 + t67 * t125;
t34 = t126 * t100 - t218;
t154 = t34 * t123 - t35 * t126;
t153 = t121 * t126 - t68 * t123;
t149 = t102 * t190 - t122 * t73;
t148 = -t117 * t207 - t120 * t124;
t147 = -g(1) * (t117 * t210 + t123 * t43) - g(2) * (-t120 * t210 + t123 * t40) - g(3) * t153;
t25 = -t120 * t211 - t126 * t40;
t27 = t117 * t211 - t126 * t43;
t146 = -g(1) * t27 - g(2) * t25 - g(3) * t51;
t145 = g(1) * t43 + g(2) * t40 - g(3) * t68;
t142 = t74 * t121;
t41 = -t117 * t152 - t120 * t142;
t44 = t117 * t142 - t120 * t152;
t144 = g(1) * t44 + g(2) * t41 - g(3) * t67;
t143 = -g(3) * t121 + (-g(1) * t117 + g(2) * t120) * t118;
t92 = pkin(2) * t180;
t141 = -pkin(2) * t212 + t41 * pkin(3) - t40 * pkin(7) + t92;
t6 = -qJDD(4) * pkin(4) - t8;
t140 = t147 - t6;
t32 = -qJD(4) * pkin(4) - t34;
t138 = -pkin(8) * t73 - t102 * t32;
t137 = -t60 * qJD(2) + t144;
t106 = -pkin(3) - t225;
t52 = -qJD(2) * pkin(3) - t54;
t136 = -qJDD(4) * t105 + (qJD(2) * t106 + t52 + t63) * qJD(4);
t135 = pkin(8) * qJD(5) * t102 + t140;
t134 = pkin(2) * t148 + t44 * pkin(3) - t43 * pkin(7);
t133 = -g(1) * t148 - g(3) * t209;
t132 = qJD(5) * t157 + t1 * t125 - t2 * t122;
t131 = -t8 * t123 + t7 * t126 + (-t123 * t35 - t126 * t34) * qJD(4);
t128 = qJD(4) ^ 2;
t28 = -qJDD(2) * pkin(3) - t30;
t130 = qJDD(2) * t106 + t105 * t128 + t137 + t28;
t91 = qJDD(4) * t126 - t128 * t123;
t90 = qJDD(4) * t123 + t128 * t126;
t81 = t159 * qJD(2);
t62 = t74 * t196;
t61 = qJD(2) * t68;
t18 = qJD(4) * t153 - t62 * t126;
t17 = qJD(4) * t51 - t62 * t123;
t14 = t122 * t81 + t125 * t34;
t13 = -t122 * t34 + t125 * t81;
t4 = qJD(5) * t19 + t61 * t122 + t18 * t125;
t3 = -qJD(5) * t20 - t18 * t122 + t61 * t125;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t200, 0, 0, 0, 0, 0, 0, (qJDD(2) * t127 - t124 * t129) * t118, (-qJDD(2) * t124 - t127 * t129) * t118, 0, -g(3) + (t121 ^ 2 + (t124 ^ 2 + t127 ^ 2) * t118 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t61 * qJD(2) - t67 * qJDD(2), t62 * qJD(2) - t68 * qJDD(2), 0, t97 * t121 - t30 * t67 + t31 * t68 - t54 * t61 - t55 * t62 - g(3), 0, 0, 0, 0, 0, 0, -t67 * t111 - t17 * qJD(4) + t153 * qJDD(4) + (-t126 * t61 + t194 * t67) * qJD(2), t67 * t184 - t18 * qJD(4) - t51 * qJDD(4) + (t123 * t61 + t193 * t67) * qJD(2), (-t123 * t153 + t126 * t51) * qJDD(2) + (t123 * t17 + t126 * t18 + (-t123 * t51 - t126 * t153) * qJD(4)) * qJD(2), t153 * t8 - t34 * t17 + t35 * t18 + t28 * t67 + t7 * t51 + t52 * t61 - g(3), 0, 0, 0, 0, 0, 0, -t3 * t102 - t153 * t47 + t17 * t76 + t19 * t73, t4 * t102 + t153 * t46 + t17 * t78 - t20 * t73, t19 * t46 - t20 * t47 - t3 * t78 - t4 * t76, t1 * t20 + t10 * t4 - t153 * t6 - t155 * t3 + t32 * t17 + t2 * t19 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t96 - g(2) * (t180 - t212) + t133, -g(1) * (t117 * t208 - t120 * t127) - g(2) * (-t117 * t127 - t120 * t208) - t200 * t124 * t118, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t119 * t213 - t137 + t30, -t116 * t213 - t145 + t201 - t31, 0, -g(2) * t92 + t54 * t60 - t55 * t63 + (g(2) * t212 + t31 * t116 + t30 * t119 + t133) * pkin(2), t114 * qJDD(2) + 0.2e1 * t161, 0.2e1 * t123 * t111 - 0.2e1 * t199 * t186, t90, t115 * qJDD(2) - 0.2e1 * t161, t91, 0, t123 * t136 - t126 * t130, t123 * t130 + t126 * t136, t131 + t145 + (qJDD(2) * t105 - t201) * (t114 + t115), -g(1) * t134 - g(2) * t141 - g(3) * t163 + t105 * t131 + t28 * t106 + t154 * t63 - t52 * t60, t78 * t173 + (-t46 * t125 - t192 * t78) * t123, -t164 + (-t181 + (t46 + t214) * t123) * t122 + t220, t167 + t227, t76 * t175 + (t123 * t47 + t193 * t76) * t122, (t47 + t177) * t126 + (t149 - t215) * t123, -t102 * t194 - t73 * t126, t48 * t73 - t221 * t102 + t145 * t122 + (-t2 + (t105 * t76 + t122 * t32) * qJD(4) - t144 * t125) * t126 + (-qJD(4) * t155 + t105 * t47 + t6 * t122 + t190 * t32 - t63 * t76) * t123, -t49 * t73 + t222 * t102 + t145 * t125 + (t1 + (t105 * t78 + t125 * t32) * qJD(4) + t144 * t122) * t126 + (-t10 * qJD(4) - t105 * t46 + t6 * t125 - t192 * t32 - t63 * t78) * t123, t48 * t46 - t49 * t47 - t221 * t78 - t222 * t76 + t157 * t193 + (-qJD(5) * t158 - t1 * t122 - t125 * t2 - t144) * t123, t1 * t49 + t2 * t48 - t32 * t123 * t63 - g(1) * (t160 * t44 + t134) - g(2) * (t160 * t41 + t141) - g(3) * (-t160 * t67 + t163) - t221 * t155 + (t6 * t123 + t193 * t32) * t105 + t222 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143 + t97, 0, 0, 0, 0, 0, 0, t91, -t90, 0, -qJD(4) * t154 + t7 * t123 + t8 * t126 + t143, 0, 0, 0, 0, 0, 0, (-t47 + t177) * t126 + (t149 + t215) * t123, t167 - t227, t164 + (t123 * t166 + t181) * t122 + t220, (qJD(4) * t158 - t6) * t126 + (qJD(4) * t32 + t132) * t123 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, t199 * t129, t184, t179, t111, qJDD(4), t88 + (-qJD(2) * t52 - t29) * t123 + t147, -t52 * t187 + (t34 + t218) * qJD(4) - t146 + t183, 0, 0, -t46 * t122 - t125 * t216, (-t46 + t217) * t125 + (-t47 + t216) * t122, (t102 * t203 - t123 * t78) * qJD(2) - t149, -t122 * t217 - t47 * t125, t102 * t192 + t125 * t73 + (-t102 * t206 + t123 * t76) * qJD(2), t102 * t195, -pkin(4) * t47 + t13 * t102 + t122 * t138 + t125 * t135 + t155 * t195 - t35 * t76, pkin(4) * t46 + t10 * t195 - t14 * t102 - t122 * t135 + t125 * t138 - t35 * t78, t13 * t78 + t14 * t76 + ((qJD(5) * t78 - t47) * pkin(8) + t232) * t125 + (pkin(8) * t166 - t228) * t122 + t146, -t10 * t14 + t155 * t13 - t32 * t35 + t140 * pkin(4) + (t132 + t146) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, -t76 ^ 2 + t78 ^ 2, -t46 - t217, -t224, -t216 - t47, t73, -t32 * t78 - g(1) * (-t27 * t122 - t44 * t125) - g(2) * (-t25 * t122 - t41 * t125) - g(3) * t19 + t228, t32 * t76 - g(1) * (t44 * t122 - t27 * t125) - g(2) * (t41 * t122 - t25 * t125) + g(3) * t20 - t232, 0, 0;];
tau_reg = t9;
