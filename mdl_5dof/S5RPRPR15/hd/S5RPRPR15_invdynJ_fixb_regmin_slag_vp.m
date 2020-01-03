% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR15
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR15_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:26
% EndTime: 2019-12-31 18:37:32
% DurationCPUTime: 2.14s
% Computational Cost: add. (1738->316), mult. (3460->433), div. (0->0), fcn. (2296->10), ass. (0->164)
t110 = sin(pkin(8));
t111 = cos(pkin(8));
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t211 = -t112 * t110 + t115 * t111;
t214 = t211 * qJD(5);
t116 = cos(qJ(3));
t177 = qJD(1) * t116;
t157 = t110 * t177;
t172 = t111 * qJD(3);
t73 = t157 - t172;
t156 = t111 * t177;
t173 = t110 * qJD(3);
t75 = t156 + t173;
t136 = t112 * t73 - t115 * t75;
t113 = sin(qJ(3));
t171 = t113 * qJD(1);
t95 = qJD(5) + t171;
t216 = t136 * t95;
t26 = t112 * t75 + t115 * t73;
t114 = sin(qJ(1));
t117 = cos(qJ(1));
t212 = g(1) * t114 - g(2) * t117;
t215 = t212 * t111;
t79 = t115 * t110 + t112 * t111;
t127 = t79 * qJD(5);
t213 = -qJD(5) + t95;
t141 = t113 * pkin(3) - t116 * qJ(4);
t208 = g(3) * t113;
t210 = t116 * t212 - t208;
t164 = t116 * qJDD(1);
t149 = -t111 * qJDD(3) + t110 * t164;
t170 = qJD(1) * qJD(3);
t155 = t113 * t170;
t42 = t110 * t155 - t149;
t199 = t110 * qJDD(3) + t111 * t164;
t43 = t111 * t155 - t199;
t7 = -t136 * qJD(5) - t112 * t43 - t115 * t42;
t118 = -pkin(1) - pkin(6);
t207 = g(3) * t116;
t205 = t26 * t95;
t154 = t116 * t170;
t165 = t113 * qJDD(1);
t128 = t154 + t165;
t77 = qJDD(5) + t128;
t204 = t211 * t77;
t203 = t79 * t77;
t202 = pkin(7) + qJ(4);
t142 = pkin(3) * t116 + qJ(4) * t113;
t58 = t142 * qJD(3) - t116 * qJD(4) + qJD(2);
t83 = qJ(2) + t141;
t20 = t58 * qJD(1) + t83 * qJDD(1);
t92 = t118 * qJD(1) + qJD(2);
t197 = t116 * t92;
t91 = t118 * qJDD(1) + qJDD(2);
t31 = qJDD(3) * qJ(4) + t113 * t91 + (qJD(4) + t197) * qJD(3);
t9 = t110 * t20 + t111 * t31;
t67 = t83 * qJD(1);
t82 = t113 * t92;
t68 = qJD(3) * qJ(4) + t82;
t23 = t110 * t67 + t111 * t68;
t130 = t211 * t113;
t201 = qJD(1) * t130 + t214;
t126 = t79 * qJD(1);
t200 = t113 * t126 + t127;
t191 = t111 * t116;
t81 = t142 * qJD(1);
t35 = t110 * t81 + t92 * t191;
t175 = qJD(3) * t116;
t158 = t118 * t175;
t30 = t110 * t58 + t111 * t158;
t189 = t113 * t118;
t41 = t110 * t83 + t111 * t189;
t176 = qJD(3) * t113;
t143 = -qJDD(3) * pkin(3) + t92 * t176 + qJDD(4);
t36 = -t116 * t91 + t143;
t196 = t36 * t110;
t195 = t36 * t111;
t194 = t36 * t116;
t193 = pkin(1) * qJDD(1);
t192 = t110 * t116;
t107 = pkin(8) + qJ(5);
t100 = sin(t107);
t188 = t114 * t100;
t101 = cos(t107);
t187 = t114 * t101;
t184 = t117 * t100;
t183 = t117 * t101;
t120 = qJD(1) ^ 2;
t182 = t120 * qJ(2);
t61 = -qJD(3) * pkin(3) + qJD(4) - t197;
t181 = -qJD(4) + t61;
t180 = t117 * pkin(1) + t114 * qJ(2);
t108 = t113 ^ 2;
t109 = t116 ^ 2;
t179 = t108 - t109;
t119 = qJD(3) ^ 2;
t178 = -t119 - t120;
t169 = qJDD(1) * qJ(2);
t168 = qJDD(1) * t110;
t167 = qJDD(1) * t111;
t166 = qJDD(3) * t113;
t163 = pkin(7) * t111 * t113;
t162 = 0.2e1 * qJD(1) * qJD(2);
t8 = -t110 * t31 + t111 * t20;
t4 = t128 * pkin(4) + t43 * pkin(7) + t8;
t5 = pkin(7) * t42 + t9;
t161 = -t112 * t5 + t115 * t4;
t160 = t110 * t171;
t159 = t113 * t173;
t153 = -t110 * t118 + pkin(4);
t152 = pkin(4) * t110 - t118;
t22 = -t110 * t68 + t111 * t67;
t34 = t111 * t81 - t92 * t192;
t148 = g(1) * t117 + g(2) * t114;
t146 = -t8 * t110 + t9 * t111;
t145 = t112 * t4 + t115 * t5;
t144 = qJDD(2) - t212;
t10 = pkin(4) * t171 - t75 * pkin(7) + t22;
t12 = -pkin(7) * t73 + t23;
t1 = t115 * t10 - t112 * t12;
t2 = t112 * t10 + t115 * t12;
t140 = t110 * t23 + t111 * t22;
t139 = -t22 * t110 + t23 * t111;
t71 = t111 * t83;
t25 = -pkin(7) * t191 + t153 * t113 + t71;
t33 = -pkin(7) * t192 + t41;
t138 = -t112 * t33 + t115 * t25;
t137 = t112 * t25 + t115 * t33;
t135 = t212 - t91;
t134 = t212 * t110;
t88 = t202 * t111;
t133 = t110 * qJD(4) + qJD(5) * t88 + (pkin(4) * t116 + t163) * qJD(1) + t34;
t87 = t202 * t110;
t132 = pkin(7) * t160 - t111 * qJD(4) + qJD(5) * t87 + t35;
t6 = -t26 * qJD(5) + t112 * t42 - t115 * t43;
t131 = t95 * t211;
t129 = 0.2e1 * qJ(2) * t170 + qJDD(3) * t118;
t125 = t135 + t182;
t123 = -t113 * t212 - t207;
t122 = -t148 + t162 + 0.2e1 * t169;
t121 = -t118 * t119 + t122;
t104 = t117 * qJ(2);
t102 = qJDD(3) * t116;
t97 = -t111 * pkin(4) - pkin(3);
t72 = t152 * t116;
t62 = t152 * t176;
t57 = t211 * t116;
t56 = t79 * t116;
t52 = t113 * t183 - t188;
t51 = t113 * t184 + t187;
t50 = t113 * t187 + t184;
t49 = -t113 * t188 + t183;
t48 = -pkin(4) * t160 + t82;
t45 = t111 * t58;
t40 = -t110 * t189 + t71;
t32 = t73 * pkin(4) + t61;
t29 = -t110 * t158 + t45;
t19 = pkin(7) * t159 + t30;
t18 = -t112 * t113 * t172 - t115 * t159 + t116 * t214;
t17 = -qJD(3) * t130 - t116 * t127;
t13 = t45 + (t153 * t116 + t163) * qJD(3);
t11 = -t42 * pkin(4) + t36;
t3 = [qJDD(1), t212, t148, t144 - 0.2e1 * t193, t122, -(qJDD(2) - t193) * pkin(1) - g(1) * (-t114 * pkin(1) + t104) - g(2) * t180 + (t162 + t169) * qJ(2), t109 * qJDD(1) - 0.2e1 * t113 * t154, -0.2e1 * t113 * t164 + 0.2e1 * t179 * t170, -t119 * t113 + t102, -t119 * t116 - t166, 0, t121 * t113 + t129 * t116, -t129 * t113 + t121 * t116, t134 + (t196 + t118 * t42 + (qJD(1) * t40 + t22) * qJD(3)) * t116 + (t29 * qJD(1) + t40 * qJDD(1) + t8 - t148 * t111 + (-t110 * t61 + t118 * t73) * qJD(3)) * t113, t215 + (t195 + t118 * t43 + (-qJD(1) * t41 - t23) * qJD(3)) * t116 + (-t30 * qJD(1) - t41 * qJDD(1) - t9 + t148 * t110 + (-t111 * t61 + t118 * t75) * qJD(3)) * t113, -t29 * t75 - t30 * t73 + t40 * t43 + t41 * t42 + t140 * t176 + (-t110 * t9 - t111 * t8 + t148) * t116, t9 * t41 + t23 * t30 + t8 * t40 + t22 * t29 - g(1) * (t117 * t141 + t104) - g(2) * (t117 * pkin(6) + t180) + (t61 * t176 - t194) * t118 + (-g(1) * t118 - g(2) * t141) * t114, -t136 * t17 + t57 * t6, t136 * t18 - t17 * t26 - t56 * t6 - t57 * t7, t6 * t113 - t136 * t175 + t17 * t95 + t57 * t77, -t7 * t113 - t175 * t26 - t18 * t95 - t56 * t77, t77 * t113 + t175 * t95, (-t112 * t19 + t115 * t13) * t95 + t138 * t77 + t161 * t113 + t1 * t175 - t62 * t26 + t72 * t7 + t11 * t56 + t32 * t18 - g(1) * t52 - g(2) * t50 + (-t113 * t2 - t137 * t95) * qJD(5), -(t112 * t13 + t115 * t19) * t95 - t137 * t77 - t145 * t113 - t2 * t175 + t62 * t136 + t72 * t6 + t11 * t57 + t32 * t17 + g(1) * t51 - g(2) * t49 + (-t1 * t113 - t138 * t95) * qJD(5); 0, 0, 0, qJDD(1), -t120, t144 - t182 - t193, 0, 0, 0, 0, 0, t178 * t113 + t102, t178 * t116 - t166, -t108 * t168 + t116 * t42 + (-t111 * t120 + (t73 - 0.2e1 * t157) * qJD(3)) * t113, -t108 * t167 + t116 * t43 + (t110 * t120 + (t75 - 0.2e1 * t156) * qJD(3)) * t113, (qJD(1) * t75 + t113 * t42 - t73 * t175) * t111 + (qJD(1) * t73 - t113 * t43 + t75 * t175) * t110, -t194 + t146 * t113 - t140 * qJD(1) + (t113 * t61 + t139 * t116) * qJD(3) - t212, 0, 0, 0, 0, 0, -qJD(1) * t131 + (-qJD(3) * t79 * t95 - t7) * t116 + (qJD(3) * t26 - t214 * t95 - t203) * t113, t95 * t126 + (-qJD(3) * t131 - t6) * t116 + (-qJD(3) * t136 + t127 * t95 - t204) * t113; 0, 0, 0, 0, 0, 0, t116 * t120 * t113, -t179 * t120, t164, -t165, qJDD(3), -t125 * t116 + t208, t125 * t113 + t207, pkin(3) * t42 - t195 + (-t215 + (-qJ(4) * t173 - t22) * qJD(1)) * t116 + (-qJ(4) * t168 + g(3) * t111 - t73 * t92 + (t181 * t110 - t34) * qJD(1)) * t113, pkin(3) * t43 + t196 + (t134 + (-qJ(4) * t172 + t23) * qJD(1)) * t116 + (-qJ(4) * t167 - g(3) * t110 - t75 * t92 + (t181 * t111 + t35) * qJD(1)) * t113, t34 * t75 + t35 * t73 + (qJ(4) * t42 - qJD(4) * t73 - t22 * t171 + t9) * t111 + (-qJ(4) * t43 + qJD(4) * t75 - t23 * t171 - t8) * t110 + t123, -t61 * t82 - t22 * t34 - t23 * t35 + t139 * qJD(4) + (-t36 - t210) * pkin(3) + (t123 + t146) * qJ(4), -t136 * t201 + t6 * t79, t136 * t200 - t201 * t26 + t211 * t6 - t7 * t79, t136 * t177 + t201 * t95 + t203, t26 * t177 - t200 * t95 + t204, -t95 * t177, (-t112 * t88 - t115 * t87) * t77 + t97 * t7 - t11 * t211 - t1 * t177 - t48 * t26 + (t112 * t132 - t115 * t133) * t95 + t200 * t32 - t210 * t101, -(-t112 * t87 + t115 * t88) * t77 + t97 * t6 + t11 * t79 + t2 * t177 + t48 * t136 + (t112 * t133 + t115 * t132) * t95 + t201 * t32 + t210 * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t75 - t173) * t171 + t149, (-t73 - t172) * t171 + t199, -t73 ^ 2 - t75 ^ 2, t135 * t116 + t22 * t75 + t23 * t73 + t143 - t208, 0, 0, 0, 0, 0, t7 - t216, t6 - t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136 * t26, t136 ^ 2 - t26 ^ 2, t6 + t205, -t7 - t216, t77, -g(1) * t49 - g(2) * t51 + t100 * t207 + t136 * t32 + t2 * t213 + t161, g(1) * t50 - g(2) * t52 + t1 * t213 + t101 * t207 + t32 * t26 - t145;];
tau_reg = t3;
