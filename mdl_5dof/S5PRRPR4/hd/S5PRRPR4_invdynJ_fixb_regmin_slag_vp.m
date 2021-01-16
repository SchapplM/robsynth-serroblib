% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:41
% EndTime: 2021-01-15 15:52:50
% DurationCPUTime: 1.94s
% Computational Cost: add. (1622->273), mult. (3709->383), div. (0->0), fcn. (2834->14), ass. (0->150)
t122 = sin(qJ(5));
t125 = cos(qJ(5));
t117 = sin(pkin(9));
t119 = cos(pkin(9));
t126 = cos(qJ(3));
t185 = t119 * t126;
t167 = qJD(2) * t185;
t123 = sin(qJ(3));
t180 = qJD(2) * t123;
t79 = t117 * t180 - t167;
t89 = t117 * t126 + t119 * t123;
t82 = t89 * qJD(2);
t152 = t122 * t79 - t125 * t82;
t171 = t126 * qJDD(2);
t172 = t123 * qJDD(2);
t156 = t117 * t172 - t119 * t171;
t81 = t89 * qJD(3);
t41 = qJD(2) * t81 + t156;
t174 = qJD(2) * qJD(3);
t166 = t123 * t174;
t137 = qJDD(2) * t89 - t117 * t166;
t165 = t126 * t174;
t42 = t119 * t165 + t137;
t132 = qJD(5) * t152 - t122 * t42 - t125 * t41;
t113 = qJD(3) + qJD(5);
t189 = t152 * t113;
t214 = t132 - t189;
t121 = qJ(4) + pkin(6);
t69 = t125 * t79;
t36 = -t122 * t82 - t69;
t188 = t36 * t113;
t178 = qJD(5) * t122;
t4 = -qJD(5) * t69 - t122 * t41 + t125 * t42 - t178 * t82;
t213 = t4 - t188;
t212 = t152 * t36;
t127 = cos(qJ(2));
t118 = sin(pkin(8));
t120 = cos(pkin(8));
t158 = g(1) * t120 + g(2) * t118;
t144 = t158 * t127;
t124 = sin(qJ(2));
t196 = g(3) * t124;
t135 = t144 + t196;
t195 = g(3) * t127;
t206 = t158 * t124;
t134 = t206 - t195;
t183 = qJDD(1) - g(3);
t211 = t127 * t183 + t206;
t210 = t152 ^ 2 - t36 ^ 2;
t114 = qJ(3) + pkin(9);
t111 = qJ(5) + t114;
t104 = sin(t111);
t105 = cos(t111);
t202 = t79 * pkin(7);
t177 = t124 * qJD(1);
t162 = qJD(2) * t121 + t177;
t75 = t162 * t126;
t190 = t119 * t75;
t191 = qJD(3) * pkin(3);
t74 = t162 * t123;
t61 = -t74 + t191;
t24 = t117 * t61 + t190;
t13 = t24 - t202;
t184 = t120 * t127;
t186 = t118 * t127;
t108 = pkin(3) * t126 + pkin(2);
t176 = t127 * qJD(1);
t85 = -qJD(2) * t108 + qJD(4) - t176;
t45 = t79 * pkin(4) + t85;
t209 = -t45 * t36 - g(1) * (-t104 * t118 - t105 * t184) - g(2) * (t104 * t120 - t105 * t186) + t13 * t178 + t105 * t196;
t175 = qJD(1) * qJD(2);
t87 = qJDD(2) * pkin(6) + t124 * qJDD(1) + t127 * t175;
t142 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + t87;
t150 = qJD(3) * t162;
t21 = qJDD(3) * pkin(3) - t123 * t142 - t126 * t150;
t22 = -t123 * t150 + t126 * t142;
t6 = -t117 * t22 + t119 * t21;
t2 = qJDD(3) * pkin(4) - pkin(7) * t42 + t6;
t7 = t117 * t21 + t119 * t22;
t3 = -pkin(7) * t41 + t7;
t208 = t45 * t152 - g(1) * (-t104 * t184 + t105 * t118) - g(2) * (-t104 * t186 - t105 * t120) - t122 * t3 + t125 * t2 + t104 * t196;
t163 = qJD(3) * t121;
t76 = t126 * qJD(4) - t123 * t163;
t77 = -t123 * qJD(4) - t126 * t163;
t194 = -t117 * t76 + t119 * t77 + t176 * t89;
t88 = t117 * t123 - t185;
t141 = t88 * t127;
t193 = qJD(1) * t141 + t117 * t77 + t119 * t76;
t84 = t88 * qJD(3);
t205 = t123 * t191 - t177;
t204 = qJD(5) - t113;
t107 = t124 * t175;
t128 = qJD(3) ^ 2;
t170 = t127 * qJDD(1);
t203 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t128 + t124 * (t158 + t175) - t107 + t170 - t195;
t201 = t82 * pkin(7);
t200 = pkin(3) * t117;
t199 = pkin(3) * t123;
t57 = t117 * t75;
t26 = -t119 * t74 - t57;
t95 = t121 * t123;
t96 = t121 * t126;
t48 = -t117 * t95 + t119 * t96;
t192 = qJD(2) * pkin(2);
t115 = t123 ^ 2;
t182 = -t126 ^ 2 + t115;
t129 = qJD(2) ^ 2;
t181 = t128 + t129;
t179 = qJD(2) * t124;
t173 = qJDD(3) * t123;
t169 = t127 * qJDD(2);
t23 = t119 * t61 - t57;
t25 = t117 * t74 - t190;
t47 = -t117 * t96 - t119 * t95;
t161 = pkin(4) * t81 + t205;
t32 = -pkin(7) * t88 + t48;
t160 = -pkin(7) * t84 + qJD(5) * t32 - t194;
t31 = -pkin(7) * t89 + t47;
t159 = pkin(7) * t81 - qJD(5) * t31 - t193;
t157 = g(1) * t118 - g(2) * t120;
t11 = qJD(3) * pkin(4) - t201 + t23;
t155 = -t122 * t11 - t125 * t13;
t72 = t89 * t124;
t73 = t88 * t124;
t154 = t122 * t73 - t125 * t72;
t153 = -t122 * t72 - t125 * t73;
t43 = t122 * t89 + t125 * t88;
t44 = -t122 * t88 + t125 * t89;
t106 = pkin(3) * t119 + pkin(4);
t149 = t106 * t122 + t125 * t200;
t148 = t106 * t125 - t122 * t200;
t143 = t157 * t126;
t136 = pkin(3) * t166 - qJDD(2) * t108 + qJDD(4) + t107;
t99 = -t176 - t192;
t131 = -pkin(6) * qJDD(3) + (t176 + t99 - t192) * qJD(3);
t46 = t136 - t170;
t130 = -qJD(2) * t99 + t135 - t87;
t112 = qJDD(3) + qJDD(5);
t110 = cos(t114);
t109 = sin(t114);
t65 = pkin(4) * t88 - t108;
t51 = pkin(3) * t180 + pkin(4) * t82;
t28 = -qJD(2) * t141 - qJD(3) * t72;
t27 = t124 * t84 - t127 * t82;
t15 = t26 - t201;
t14 = t25 + t202;
t12 = t41 * pkin(4) + t46;
t9 = qJD(5) * t44 - t122 * t84 + t125 * t81;
t8 = -qJD(5) * t43 - t122 * t81 - t125 * t84;
t1 = [t183, 0, -t124 * t129 + t169, -qJDD(2) * t124 - t127 * t129, 0, 0, 0, 0, 0, (-0.2e1 * t166 + t171) * t127 + (-t126 * t181 - t173) * t124, (-qJDD(3) * t124 - 0.2e1 * t127 * t174) * t126 + (t124 * t181 - t169) * t123, qJD(3) * t27 - qJDD(3) * t72 - t127 * t41 + t179 * t79, -qJD(3) * t28 + qJDD(3) * t73 - t127 * t42 + t179 * t82, -t27 * t82 - t28 * t79 + t41 * t73 + t42 * t72, -t127 * t46 + t179 * t85 + t23 * t27 + t24 * t28 - t6 * t72 - t7 * t73 - g(3), 0, 0, 0, 0, 0, (-qJD(5) * t153 - t122 * t28 + t125 * t27) * t113 + t154 * t112 - t36 * t179 + t127 * t132, -(qJD(5) * t154 + t122 * t27 + t125 * t28) * t113 - t153 * t112 - t152 * t179 - t127 * t4; 0, qJDD(2), t211, -t124 * t183 + t144, qJDD(2) * t115 + 0.2e1 * t123 * t165, 0.2e1 * t123 * t171 - 0.2e1 * t174 * t182, t126 * t128 + t173, qJDD(3) * t126 - t123 * t128, 0, t123 * t131 + t126 * t203, -t123 * t203 + t126 * t131, -t79 * t177 + t47 * qJDD(3) - t108 * t41 + t46 * t88 + t85 * t81 + t134 * t110 + (t199 * t79 + t194) * qJD(3), -t82 * t177 - t48 * qJDD(3) - t108 * t42 + t46 * t89 - t85 * t84 - t134 * t109 + (t199 * t82 - t193) * qJD(3), -t193 * t79 - t194 * t82 + t23 * t84 - t24 * t81 - t48 * t41 - t47 * t42 - t6 * t89 - t7 * t88 - t135, t7 * t48 + t6 * t47 - t46 * t108 - g(3) * (t108 * t127 + t121 * t124) + t205 * t85 + t193 * t24 + t194 * t23 + t158 * (t108 * t124 - t121 * t127), -t152 * t8 + t4 * t44, t132 * t44 + t152 * t9 + t36 * t8 - t4 * t43, t112 * t44 + t113 * t8, -t112 * t43 - t113 * t9, 0, (-t122 * t32 + t125 * t31) * t112 - t65 * t132 + t12 * t43 + t45 * t9 - t161 * t36 + (t122 * t159 - t125 * t160) * t113 + t134 * t105, -(t122 * t31 + t125 * t32) * t112 + t65 * t4 + t12 * t44 + t45 * t8 - t161 * t152 + (t122 * t160 + t125 * t159) * t113 - t134 * t104; 0, 0, 0, 0, -t123 * t129 * t126, t182 * t129, t172, t171, qJDD(3), t123 * t130 - t143, t123 * t157 + t126 * t130, -t25 * qJD(3) - t85 * t82 - g(1) * (-t109 * t184 + t110 * t118) - g(2) * (-t109 * t186 - t110 * t120) + t109 * t196 + (qJDD(3) * t119 - t180 * t79) * pkin(3) + t6, t26 * qJD(3) + t85 * t79 - g(1) * (-t109 * t118 - t110 * t184) - g(2) * (t109 * t120 - t110 * t186) + t110 * t196 + (-qJDD(3) * t117 - t180 * t82) * pkin(3) - t7, (t24 + t25) * t82 + (-t23 + t26) * t79 + (-t117 * t41 - t119 * t42) * pkin(3), -t23 * t25 - t24 * t26 + (t7 * t117 + t6 * t119 - t143 + (-t85 * qJD(2) + t135) * t123) * pkin(3), t212, t210, t213, t214, t112, t148 * t112 - (-t122 * t15 + t125 * t14) * t113 + t51 * t36 + (-t113 * t149 + t155) * qJD(5) + t208, -t149 * t112 - t125 * t3 - t122 * t2 + (t122 * t14 + t125 * t15) * t113 + t51 * t152 + (-t125 * t11 - t113 * t148) * qJD(5) + t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t82 * qJD(3) + t156, (-t79 + t167) * qJD(3) + t137, -t79 ^ 2 - t82 ^ 2, t23 * t82 + t24 * t79 + t136 - t211, 0, 0, 0, 0, 0, -t132 - t189, t4 + t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t210, t213, t214, t112, t155 * t204 + t208, (-t113 * t13 - t2) * t122 + (-t11 * t204 - t3) * t125 + t209;];
tau_reg = t1;
