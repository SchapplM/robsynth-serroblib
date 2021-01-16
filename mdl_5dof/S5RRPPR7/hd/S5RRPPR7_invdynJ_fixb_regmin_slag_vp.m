% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:00
% EndTime: 2021-01-15 19:59:09
% DurationCPUTime: 2.22s
% Computational Cost: add. (1989->316), mult. (4644->397), div. (0->0), fcn. (3265->10), ass. (0->164)
t119 = cos(qJ(5));
t113 = sin(pkin(8));
t114 = cos(pkin(8));
t117 = sin(qJ(2));
t120 = cos(qJ(2));
t82 = t113 * t120 + t114 * t117;
t207 = t82 * qJD(1);
t210 = qJD(5) + t207;
t116 = sin(qJ(5));
t211 = t116 * t210;
t163 = qJD(1) * qJD(2);
t156 = t117 * t163;
t132 = qJDD(1) * t82 - t113 * t156;
t155 = t120 * t163;
t46 = t114 * t155 + t132;
t41 = qJDD(5) + t46;
t139 = t119 * t41 - t210 * t211;
t166 = t116 * qJD(2);
t180 = t114 * t120;
t157 = qJD(1) * t180;
t169 = qJD(1) * t117;
t71 = t113 * t169 - t157;
t54 = -t119 * t71 + t166;
t212 = t210 * t54;
t161 = t120 * qJDD(1);
t162 = t117 * qJDD(1);
t141 = t113 * t162 - t114 * t161;
t209 = 0.2e1 * t207 * qJD(2) + t141;
t70 = t207 ^ 2;
t208 = -t71 ^ 2 - t70;
t115 = -qJ(3) - pkin(6);
t90 = t115 * t120;
t85 = qJD(1) * t90;
t77 = t113 * t85;
t89 = t115 * t117;
t84 = qJD(1) * t89;
t49 = t114 * t84 + t77;
t173 = -qJD(4) + t49;
t206 = -qJD(5) + t210;
t118 = sin(qJ(1));
t121 = cos(qJ(1));
t205 = g(1) * t118 - g(2) * t121;
t197 = t71 * pkin(4);
t186 = t114 * t85;
t80 = qJD(2) * pkin(2) + t84;
t44 = t113 * t80 - t186;
t35 = -qJD(2) * qJ(4) - t44;
t21 = -t35 - t197;
t48 = t113 * t84 - t186;
t100 = -t114 * pkin(2) - pkin(3);
t95 = -pkin(7) + t100;
t204 = t95 * t41 + (t21 - t48 + t197) * t210;
t109 = qJ(2) + pkin(8);
t105 = sin(t109);
t53 = t113 * t89 - t114 * t90;
t203 = -t53 * qJDD(2) - t105 * t205;
t145 = g(1) * t121 + g(2) * t118;
t202 = -t48 * qJD(2) - t145 * t105;
t106 = cos(t109);
t52 = -t113 * t90 - t114 * t89;
t201 = -t52 * qJDD(2) + t106 * t205;
t130 = -g(3) * t105 - t145 * t106;
t199 = pkin(3) + pkin(7);
t73 = t82 * qJD(2);
t45 = qJD(1) * t73 + t141;
t198 = t45 * pkin(3);
t196 = t207 * pkin(4);
t195 = pkin(2) * t117;
t98 = g(3) * t106;
t192 = g(3) * t120;
t191 = t120 * pkin(2);
t102 = pkin(1) + t191;
t140 = -t82 * qJ(4) - t102;
t81 = t113 * t117 - t180;
t22 = t199 * t81 + t140;
t190 = t22 * t41;
t56 = t119 * qJD(2) + t116 * t71;
t189 = t56 * t71;
t188 = t71 * t54;
t151 = qJD(2) * t115;
t69 = -t117 * qJD(3) + t120 * t151;
t40 = qJDD(2) * pkin(2) + qJD(1) * t69 + qJDD(1) * t89;
t68 = t120 * qJD(3) + t117 * t151;
t47 = qJD(1) * t68 - qJDD(1) * t90;
t187 = t113 * t47 - t114 * t40;
t11 = t113 * t40 + t114 * t47;
t167 = qJD(5) * t119;
t13 = -qJD(5) * t166 + t119 * qJDD(2) + t116 * t45 + t71 * t167;
t184 = t13 * t119;
t183 = t46 * qJ(4);
t182 = qJD(5) * t81;
t181 = qJDD(2) * pkin(3);
t179 = t118 * t115;
t178 = t118 * t116;
t177 = t118 * t119;
t176 = t121 * t116;
t175 = t121 * t119;
t172 = t196 - t173;
t111 = t117 ^ 2;
t170 = -t120 ^ 2 + t111;
t168 = qJD(2) * t117;
t160 = qJDD(2) * qJ(4) + t11;
t104 = pkin(2) * t168;
t159 = t81 * t167;
t158 = qJDD(4) + t187;
t61 = pkin(2) * t156 - qJDD(1) * t102 + qJDD(3);
t126 = -qJD(4) * t207 - t183 + t61;
t1 = t199 * t45 + t126;
t43 = t114 * t80 + t77;
t146 = qJD(4) - t43;
t16 = -t199 * qJD(2) + t146 + t196;
t153 = qJD(5) * t16 + t1;
t87 = -qJD(1) * t102 + qJD(3);
t131 = -qJ(4) * t207 + t87;
t15 = t199 * t71 + t131;
t5 = t46 * pkin(4) - t199 * qJDD(2) + t158;
t152 = -qJD(5) * t15 + t5;
t27 = t113 * t68 - t114 * t69;
t150 = pkin(2) * t169 + t71 * qJ(4);
t148 = t116 * qJDD(2) - t119 * t45;
t142 = t210 * t73 + t41 * t81;
t28 = t113 * t69 + t114 * t68;
t3 = t116 * t16 + t119 * t15;
t8 = -qJD(2) * qJD(4) - t160;
t76 = qJD(2) * t180 - t113 * t168;
t138 = -t76 * qJ(4) - t82 * qJD(4) + t104;
t137 = -0.2e1 * pkin(1) * t163 - pkin(6) * qJDD(2);
t25 = t71 * pkin(3) + t131;
t136 = t207 * t25 + t158 + t98;
t29 = t82 * pkin(4) + t52;
t6 = -t45 * pkin(4) - t8;
t135 = t21 * t73 - t29 * t41 + t6 * t81;
t133 = -t119 * t210 ^ 2 - t116 * t41;
t122 = qJD(2) ^ 2;
t129 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t122 + t205;
t123 = qJD(1) ^ 2;
t128 = pkin(1) * t123 - pkin(6) * qJDD(1) + t145;
t127 = -t205 + t61;
t125 = t207 * t27 - t28 * t71 - t53 * t45 + t52 * t46 - t145;
t124 = t6 + (-qJD(5) * t95 + t199 * t207 + t150) * t210 + t130;
t101 = t115 * t121;
t97 = t113 * pkin(2) + qJ(4);
t88 = -t113 * pkin(3) + qJ(4) * t114;
t86 = pkin(3) * t114 + qJ(4) * t113 + pkin(2);
t67 = -t105 * t178 + t175;
t66 = t105 * t177 + t176;
t65 = t105 * t176 + t177;
t64 = t105 * t175 - t178;
t62 = qJD(2) * t71;
t50 = t88 * t117 + t86 * t120 + pkin(1);
t42 = t81 * pkin(3) + t140;
t31 = -qJD(2) * pkin(3) + t146;
t30 = -t81 * pkin(4) + t53;
t26 = pkin(3) * t207 + t150;
t20 = t73 * pkin(3) + t138;
t19 = -t73 * pkin(4) + t28;
t18 = t76 * pkin(4) + t27;
t14 = qJD(5) * t56 + t148;
t12 = t199 * t73 + t138;
t9 = t158 - t181;
t7 = t126 + t198;
t4 = t119 * t5;
t2 = -t116 * t15 + t119 * t16;
t10 = [qJDD(1), t205, t145, t111 * qJDD(1) + 0.2e1 * t117 * t155, 0.2e1 * t117 * t161 - 0.2e1 * t163 * t170, qJDD(2) * t117 + t122 * t120, qJDD(2) * t120 - t122 * t117, 0, t117 * t137 + t120 * t129, -t117 * t129 + t120 * t137, -t102 * t45 + t61 * t81 + t87 * t73 + (t71 * t195 - t27) * qJD(2) + t201, -t102 * t46 + t61 * t82 + t87 * t76 + (t195 * t207 - t28) * qJD(2) + t203, -t11 * t81 + t187 * t82 - t43 * t76 - t44 * t73 + t125, t11 * t53 + t44 * t28 + t187 * t52 - t43 * t27 - t61 * t102 + t87 * t104 - g(1) * (-t118 * t102 - t101) - g(2) * (t121 * t102 - t179), t31 * t76 + t35 * t73 + t8 * t81 + t9 * t82 + t125, t27 * qJD(2) - t20 * t71 - t25 * t73 - t42 * t45 - t7 * t81 - t201, t28 * qJD(2) - t20 * t207 - t25 * t76 - t42 * t46 - t7 * t82 - t203, t7 * t42 + t25 * t20 - t8 * t53 - t35 * t28 + t9 * t52 + t31 * t27 - g(1) * (-t50 * t118 - t101) - g(2) * (t50 * t121 - t179), t56 * t159 + (t13 * t81 + t56 * t73) * t116, (-t116 * t54 + t119 * t56) * t73 + (-t116 * t14 + t184 + (-t116 * t56 - t119 * t54) * qJD(5)) * t81, t116 * t142 + t13 * t82 + t159 * t210 + t56 * t76, t119 * t142 - t14 * t82 - t182 * t211 - t54 * t76, t210 * t76 + t41 * t82, -g(1) * t67 - g(2) * t65 + t30 * t14 + t19 * t54 + t2 * t76 + t4 * t82 + (-t1 * t82 - t12 * t210 - t190) * t116 + (t18 * t210 - t135) * t119 + ((-t116 * t29 - t119 * t22) * t210 - t3 * t82 + t21 * t116 * t81) * qJD(5), g(1) * t66 - g(2) * t64 + t30 * t13 + t19 * t56 - t3 * t76 + (-(qJD(5) * t29 + t12) * t210 - t190 - t153 * t82 + t21 * t182) * t119 + (-(-qJD(5) * t22 + t18) * t210 - t152 * t82 + t135) * t116; 0, 0, 0, -t117 * t123 * t120, t170 * t123, t162, t161, qJDD(2), t117 * t128 - t192, g(3) * t117 + t120 * t128, -t87 * t207 - t98 + (qJDD(2) * t114 - t169 * t71) * pkin(2) - t187 - t202, t49 * qJD(2) + t87 * t71 + (-qJDD(2) * t113 - t169 * t207) * pkin(2) - t11 - t130, (t44 - t48) * t207 + (-t43 + t49) * t71 + (-t113 * t45 - t114 * t46) * pkin(2), t43 * t48 - t44 * t49 + (-t192 - t187 * t114 + t11 * t113 + (-qJD(1) * t87 + t145) * t117) * pkin(2), t100 * t46 - t97 * t45 + (-t35 - t48) * t207 + (t31 + t173) * t71, t26 * t71 + (-pkin(3) + t100) * qJDD(2) + t136 + t202, t97 * qJDD(2) - t25 * t71 + t26 * t207 + (0.2e1 * qJD(4) - t49) * qJD(2) + t130 + t160, -t8 * t97 + t9 * t100 - t25 * t26 - t31 * t48 - g(3) * (t106 * pkin(3) + t105 * qJ(4) + t191) - t145 * (-t86 * t117 + t88 * t120) + t173 * t35, -t211 * t56 + t184, (-t210 * t56 - t14) * t119 + (-t13 + t212) * t116, t139 + t189, t133 - t188, t210 * t71, t124 * t116 + t204 * t119 + t97 * t14 + t172 * t54 + t2 * t71, -t204 * t116 + t124 * t119 + t97 * t13 + t172 * t56 - t3 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, (-t71 + t157) * qJD(2) + t132, t208, t207 * t43 + t44 * t71 + t127, t208, -t209, t62 - t46, t198 - t183 - t35 * t71 + (-qJD(4) - t31) * t207 + t127, 0, 0, 0, 0, 0, t133 + t188, -t139 + t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 + t46, -t207 * t71 + qJDD(2), -t70 - t122, t35 * qJD(2) - t145 * t82 + t136 - t181, 0, 0, 0, 0, 0, -qJD(2) * t54 + t139, -qJD(2) * t56 + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, t13 + t212, t206 * t56 - t148, t41, -g(1) * t64 - g(2) * t66 - t116 * t1 + t119 * t98 + t206 * t3 - t21 * t56 + t4, g(1) * t65 - g(2) * t67 + t2 * t210 + t21 * t54 - t153 * t119 + (-t152 - t98) * t116;];
tau_reg = t10;
