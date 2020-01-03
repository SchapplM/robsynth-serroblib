% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:58
% EndTime: 2019-12-31 18:28:01
% DurationCPUTime: 1.52s
% Computational Cost: add. (1568->255), mult. (3652->314), div. (0->0), fcn. (2791->10), ass. (0->150)
t113 = qJD(3) - qJD(5);
t120 = sin(qJ(5));
t123 = cos(qJ(5));
t118 = cos(pkin(8));
t203 = cos(qJ(3));
t172 = t203 * t118;
t159 = qJD(1) * t172;
t117 = sin(pkin(8));
t121 = sin(qJ(3));
t184 = t121 * t117;
t171 = qJD(1) * t184;
t75 = -t159 + t171;
t85 = t203 * t117 + t121 * t118;
t77 = t85 * qJD(1);
t211 = t120 * t75 + t123 * t77;
t193 = t211 * t113;
t166 = qJDD(1) * t203;
t175 = t118 * qJDD(1);
t173 = qJD(3) * t159 + t117 * t166 + t121 * t175;
t41 = qJD(3) * t171 - t173;
t176 = t117 * qJDD(1);
t155 = -t118 * t166 + t121 * t176;
t80 = t85 * qJD(3);
t42 = qJD(1) * t80 + t155;
t5 = qJD(5) * t211 - t120 * t41 - t123 * t42;
t220 = t5 + t193;
t104 = t118 * pkin(2) + pkin(1);
t88 = -t104 * qJDD(1) + qJDD(2);
t219 = t41 * qJ(4) + t88;
t89 = -t104 * qJD(1) + qJD(2);
t218 = -t77 * qJ(4) + t89;
t198 = pkin(6) + qJ(2);
t92 = t198 * t118;
t87 = qJD(1) * t92;
t72 = t121 * t87;
t91 = t198 * t117;
t86 = qJD(1) * t91;
t45 = -t203 * t86 - t72;
t180 = qJD(4) - t45;
t188 = qJDD(1) * pkin(1);
t124 = cos(qJ(1));
t122 = sin(qJ(1));
t202 = g(1) * t122;
t212 = -g(2) * t124 + t202;
t139 = -qJDD(2) + t188 + t212;
t148 = t120 * t77 - t123 * t75;
t164 = t148 * qJD(5) - t120 * t42 + t123 * t41;
t192 = t148 * t113;
t217 = t164 + t192;
t216 = -t148 ^ 2 + t211 ^ 2;
t125 = -pkin(3) - pkin(4);
t12 = t125 * t75 - t218;
t112 = pkin(8) + qJ(3);
t106 = sin(t112);
t107 = cos(t112);
t142 = t106 * t120 + t107 * t123;
t170 = qJD(3) * t203;
t178 = qJD(3) * t121;
t177 = qJD(1) * qJD(2);
t206 = t198 * qJDD(1) + t177;
t59 = t206 * t117;
t60 = t206 * t118;
t165 = t121 * t60 + t87 * t170 - t86 * t178 + t203 * t59;
t156 = qJDD(4) + t165;
t2 = t41 * pkin(7) + t125 * qJDD(3) + t156;
t114 = qJDD(3) * qJ(4);
t115 = qJD(3) * qJD(4);
t174 = -t121 * t59 - t86 * t170 + t203 * t60;
t7 = -t87 * t178 + t114 + t115 + t174;
t3 = t42 * pkin(7) + t7;
t68 = t106 * t123 - t107 * t120;
t53 = t68 * t122;
t185 = t107 * t124;
t186 = t106 * t124;
t55 = t120 * t185 - t123 * t186;
t215 = -g(1) * t55 + g(2) * t53 - g(3) * t142 + t12 * t211 + t120 * t3 - t123 * t2;
t213 = t211 * t148;
t200 = g(2) * t122;
t158 = g(1) * t124 + t200;
t181 = -t77 * pkin(7) + t180;
t210 = qJD(5) + t113;
t209 = qJ(2) * qJDD(1);
t54 = t142 * t122;
t56 = t142 * t124;
t208 = g(1) * t56 + g(2) * t54 + g(3) * t68 + t12 * t148 - t120 * t2 - t123 * t3;
t25 = qJD(2) * t172 - t91 * t170 + (-qJD(2) * t117 - qJD(3) * t92) * t121;
t48 = -t121 * t91 + t203 * t92;
t207 = -t25 * qJD(3) - t48 * qJDD(3) - t106 * t212;
t205 = g(3) * t106 + qJD(3) * (t45 + t72) + t107 * t158 - t174;
t204 = t77 ^ 2;
t199 = t77 * t75;
t46 = -t121 * t86 + t203 * t87;
t190 = t75 * qJ(4);
t187 = qJDD(3) * pkin(3);
t183 = t46 * qJD(3);
t182 = t77 * qJD(4);
t179 = t117 ^ 2 + t118 ^ 2;
t163 = t179 * qJD(1) ^ 2;
t162 = t113 ^ 2;
t23 = t75 * pkin(7) + t46;
t160 = 0.2e1 * t179;
t154 = t107 * pkin(3) + t106 * qJ(4);
t15 = t125 * qJD(3) + t181;
t116 = qJD(3) * qJ(4);
t16 = t116 + t23;
t152 = t120 * t16 - t123 * t15;
t151 = -t120 * t15 - t123 * t16;
t47 = t121 * t92 + t203 * t91;
t27 = -t85 * pkin(7) + t47;
t84 = -t172 + t184;
t28 = t84 * pkin(7) + t48;
t150 = -t120 * t28 + t123 * t27;
t149 = t120 * t27 + t123 * t28;
t43 = t120 * t85 - t123 * t84;
t44 = t120 * t84 + t123 * t85;
t79 = t117 * t178 - t118 * t170;
t145 = -t79 * qJ(4) + t85 * qJD(4);
t144 = t123 * qJ(4) + t120 * t125;
t143 = -t120 * qJ(4) + t123 * t125;
t140 = t85 * qJ(4) + t104;
t138 = t104 + t154;
t137 = g(1) * t186 - g(3) * t107 + t106 * t200 - t165;
t136 = t139 + t188;
t26 = t85 * qJD(2) + t48 * qJD(3);
t135 = -g(2) * t185 - t26 * qJD(3) - t47 * qJDD(3) + t107 * t202;
t24 = t75 * pkin(3) + t218;
t131 = t24 * t77 + qJDD(4) - t137;
t130 = t42 * pkin(3) + t219;
t129 = t160 * t177 - t158;
t128 = 0.2e1 * t77 * qJD(3) + t155;
t109 = qJDD(3) - qJDD(5);
t69 = t75 ^ 2;
t40 = t84 * pkin(3) - t140;
t39 = t77 * pkin(3) + t190;
t38 = t116 + t46;
t37 = -qJD(3) * pkin(3) + t180;
t21 = t125 * t84 + t140;
t20 = t80 * pkin(3) - t145;
t19 = (t75 - t171) * qJD(3) + t173;
t18 = (t75 + t171) * qJD(3) - t173;
t17 = t125 * t77 - t190;
t14 = t79 * pkin(7) + t26;
t13 = t80 * pkin(7) + t25;
t11 = t125 * t80 + t145;
t10 = t44 * qJD(5) - t120 * t79 - t123 * t80;
t9 = -t43 * qJD(5) + t120 * t80 - t123 * t79;
t8 = t156 - t187;
t6 = t130 - t182;
t1 = t125 * t42 + t182 - t219;
t4 = [qJDD(1), t212, t158, t136 * t118, -t136 * t117, t160 * t209 + t129, pkin(1) * t139 + (t179 * t209 + t129) * qJ(2), -t41 * t85 - t77 * t79, t41 * t84 - t85 * t42 + t79 * t75 - t77 * t80, -t79 * qJD(3) + t85 * qJDD(3), -t80 * qJD(3) - t84 * qJDD(3), 0, -t104 * t42 + t89 * t80 + t88 * t84 + t135, t104 * t41 - t89 * t79 + t88 * t85 + t207, t20 * t75 + t24 * t80 + t40 * t42 + t6 * t84 + t135, -t25 * t75 + t26 * t77 - t37 * t79 - t38 * t80 - t47 * t41 - t48 * t42 - t7 * t84 + t8 * t85 - t158, -t20 * t77 + t24 * t79 + t40 * t41 - t6 * t85 - t207, t24 * t20 + t38 * t25 + t37 * t26 + t6 * t40 + t8 * t47 + t7 * t48 + (-g(1) * t198 - g(2) * t138) * t124 + (g(1) * t138 - g(2) * t198) * t122, -t164 * t44 + t211 * t9, -t10 * t211 - t148 * t9 + t164 * t43 - t44 * t5, -t44 * t109 - t9 * t113, t10 * t113 + t43 * t109, 0, t11 * t148 + t21 * t5 + t1 * t43 + t12 * t10 - (-qJD(5) * t149 - t120 * t13 + t123 * t14) * t113 - t150 * t109 + g(1) * t54 - g(2) * t56, t11 * t211 - t21 * t164 + t1 * t44 + t12 * t9 + (qJD(5) * t150 + t120 * t14 + t123 * t13) * t113 + t149 * t109 + g(1) * t53 + g(2) * t55; 0, 0, 0, -t175, t176, -t163, -qJ(2) * t163 - t139, 0, 0, 0, 0, 0, t128, -t18, t128, -t69 - t204, t18, t38 * t75 + (-qJD(4) - t37) * t77 + t130 - t212, 0, 0, 0, 0, 0, -t5 + t193, t164 - t192; 0, 0, 0, 0, 0, 0, 0, t199, -t69 + t204, t19, -t155, qJDD(3), -t89 * t77 + t137 + t183, t89 * t75 + t205, -t39 * t75 - t131 + t183 + 0.2e1 * t187, pkin(3) * t41 - t42 * qJ(4) + (t38 - t46) * t77 + (t37 - t180) * t75, -t24 * t75 + t39 * t77 + 0.2e1 * t114 + 0.2e1 * t115 - t205, -t8 * pkin(3) - g(3) * t154 + t7 * qJ(4) + t180 * t38 - t24 * t39 - t37 * t46 + t158 * (pkin(3) * t106 - qJ(4) * t107), -t213, -t216, t217, t220, t109, -t143 * t109 - t17 * t148 + (t120 * t181 + t123 * t23) * t113 + (t113 * t144 - t151) * qJD(5) + t215, t144 * t109 - t17 * t211 + (-t120 * t23 + t123 * t181) * t113 + (t113 * t143 - t152) * qJD(5) - t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t199, t19, -qJD(3) ^ 2 - t204, -t38 * qJD(3) + t131 - t187, 0, 0, 0, 0, 0, -t123 * t109 - t120 * t162 - t148 * t77, t120 * t109 - t123 * t162 - t211 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t216, -t217, -t220, -t109, t151 * t210 - t215, t152 * t210 + t208;];
tau_reg = t4;
