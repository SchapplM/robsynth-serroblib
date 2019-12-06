% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:29
% EndTime: 2019-12-05 16:56:36
% DurationCPUTime: 2.10s
% Computational Cost: add. (1659->300), mult. (3945->436), div. (0->0), fcn. (3004->10), ass. (0->159)
t104 = sin(qJ(3));
t166 = qJD(2) * t104;
t218 = qJD(4) * t166 - qJDD(3);
t103 = sin(qJ(4));
t106 = cos(qJ(4));
t154 = t104 * qJDD(2);
t107 = cos(qJ(3));
t165 = qJD(2) * t107;
t23 = ((qJD(4) + t165) * qJD(3) + t154) * t103 + t218 * t106;
t92 = pkin(4) * t106 + pkin(3);
t217 = -t107 * t92 - pkin(2);
t151 = pkin(4) * t103 + pkin(7);
t105 = sin(qJ(2));
t163 = qJD(3) * t104;
t100 = sin(pkin(5));
t169 = qJD(1) * t100;
t108 = cos(qJ(2));
t171 = t107 * t108;
t201 = pkin(7) * t103;
t133 = pkin(3) * t104 - pkin(8) * t107;
t72 = t133 * qJD(3);
t216 = (-t103 * t171 + t105 * t106) * t169 - t106 * t72 - t163 * t201;
t160 = qJD(4) * t106;
t75 = -pkin(3) * t107 - pkin(8) * t104 - pkin(2);
t215 = -(t103 * t105 + t106 * t171) * t169 + t103 * t72 + t75 * t160;
t101 = cos(pkin(5));
t174 = t101 * t107;
t73 = qJD(2) * pkin(7) + t105 * t169;
t214 = qJD(1) * t174 - t104 * t73;
t183 = cos(pkin(9));
t139 = t100 * t183;
t138 = t183 * t105;
t99 = sin(pkin(9));
t184 = t99 * t108;
t50 = t101 * t138 + t184;
t26 = -t104 * t139 + t107 * t50;
t137 = t183 * t108;
t188 = t105 * t99;
t52 = -t101 * t188 + t137;
t28 = t100 * t104 * t99 + t107 * t52;
t175 = t100 * t108;
t176 = t100 * t107;
t55 = t101 * t104 + t105 * t176;
t31 = -t103 * t55 - t106 * t175;
t49 = -t101 * t137 + t188;
t51 = t101 * t184 + t138;
t213 = -g(1) * (-t103 * t28 + t106 * t51) - g(2) * (-t103 * t26 + t106 * t49) - g(3) * t31;
t212 = pkin(4) * t23 + qJDD(5);
t109 = qJD(3) ^ 2;
t157 = qJD(1) * qJD(2);
t145 = t105 * t157;
t132 = -qJDD(1) * t175 + t100 * t145;
t135 = g(1) * t51 + g(2) * t49;
t211 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t109 + t100 * (-g(3) * t108 + t145) - t132 + t135;
t118 = g(3) * t175 - t135;
t168 = qJD(1) * t104;
t86 = t101 * t168;
t43 = t107 * t73 + t86;
t37 = qJD(3) * pkin(8) + t43;
t88 = -qJD(4) + t165;
t209 = (pkin(7) * t88 + t37) * qJD(4) - t118;
t164 = qJD(3) * t103;
t69 = t106 * t166 + t164;
t208 = t69 ^ 2;
t149 = t108 * t169;
t44 = t75 * qJD(2) - t149;
t15 = -t103 * t37 + t106 * t44;
t9 = -qJ(5) * t69 + t15;
t4 = -pkin(4) * t88 + t9;
t207 = -t9 + t4;
t172 = t106 * t107;
t129 = pkin(4) * t104 - qJ(5) * t172;
t158 = t106 * qJD(5);
t182 = qJ(5) * t104;
t89 = pkin(7) * t172;
t204 = -t104 * t158 + t129 * qJD(3) + (-t89 + (-t75 + t182) * t103) * qJD(4) - t216;
t173 = t104 * t106;
t203 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t173 + (-qJD(5) * t104 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t107) * t103 + t215;
t159 = t106 * qJD(3);
t67 = t103 * t166 - t159;
t200 = t67 * t88;
t199 = t69 * t88;
t198 = qJ(5) + pkin(8);
t71 = t133 * qJD(2);
t197 = t103 * t71 + t106 * t214;
t140 = qJD(4) * t198;
t196 = t158 - t197 + (qJ(5) * t165 - t140) * t103;
t58 = t106 * t71;
t195 = -t129 * qJD(2) - t106 * t140 - t58 + (-qJD(5) + t214) * t103;
t192 = t103 * t75 + t89;
t97 = t104 ^ 2;
t191 = -t107 ^ 2 + t97;
t190 = qJD(2) * pkin(2);
t189 = t103 * t88;
t187 = t106 * t69;
t156 = qJD(2) * qJD(3);
t144 = t107 * t156;
t22 = -qJD(4) * t159 + (-t144 - t154) * t106 + t218 * t103;
t185 = t22 * t103;
t181 = qJD(3) * t67;
t180 = qJD(4) * t88;
t178 = qJDD(3) * pkin(3);
t177 = t100 * t105;
t170 = qJDD(1) - g(3);
t167 = qJD(2) * t100;
t162 = qJD(3) * t107;
t161 = qJD(4) * t103;
t155 = qJDD(1) * t101;
t95 = t107 * qJDD(2);
t141 = t104 * t155;
t46 = qJDD(2) * pkin(7) + (qJDD(1) * t105 + t108 * t157) * t100;
t12 = qJDD(3) * pkin(8) + qJD(3) * t214 + t107 * t46 + t141;
t20 = qJD(2) * t72 + t75 * qJDD(2) + t132;
t153 = t103 * t20 + t106 * t12 + t44 * t160;
t152 = t88 * t159;
t147 = t105 * t167;
t146 = t108 * t167;
t143 = t108 * t156;
t134 = g(1) * t52 + g(2) * t50;
t16 = t103 * t44 + t106 * t37;
t110 = qJD(2) ^ 2;
t131 = qJDD(2) * t108 - t105 * t110;
t64 = t104 * t156 + qJDD(4) - t95;
t128 = t103 * t64 - t88 * t160;
t127 = t106 * t64 + t88 * t161;
t126 = t103 * t175 - t106 * t55;
t54 = t104 * t177 - t174;
t123 = t37 * t161 - t153;
t25 = t50 * t104 + t107 * t139;
t27 = t104 * t52 - t99 * t176;
t122 = g(1) * t27 + g(2) * t25 + g(3) * t54;
t121 = g(1) * t28 + g(2) * t26 + g(3) * t55;
t120 = qJD(3) * t86 + t104 * t46 - t107 * t155 + t73 * t162;
t36 = -qJD(3) * pkin(3) - t214;
t119 = -g(3) * t177 - t134;
t117 = -pkin(8) * t64 - t88 * t36;
t13 = t120 - t178;
t19 = t106 * t20;
t115 = -t16 * qJD(4) - t103 * t12 + t19;
t114 = pkin(8) * t180 + t122 - t13;
t74 = -t149 - t190;
t113 = -pkin(7) * qJDD(3) + (t149 + t74 - t190) * qJD(3);
t112 = -t120 + t122;
t77 = t198 * t106;
t76 = t198 * t103;
t66 = t106 * t75;
t63 = t67 ^ 2;
t33 = -t103 * t182 + t192;
t30 = t55 * qJD(3) + t104 * t146;
t29 = -t54 * qJD(3) + t107 * t146;
t24 = -qJ(5) * t173 + t66 + (-pkin(4) - t201) * t107;
t21 = pkin(4) * t67 + qJD(5) + t36;
t10 = -qJ(5) * t67 + t16;
t7 = t31 * qJD(4) + t103 * t147 + t106 * t29;
t6 = t126 * qJD(4) - t103 * t29 + t106 * t147;
t3 = t13 + t212;
t2 = -qJ(5) * t23 - qJD(5) * t67 - t123;
t1 = pkin(4) * t64 + qJ(5) * t22 - qJD(5) * t69 + t115;
t5 = [t170, 0, t131 * t100, (-qJDD(2) * t105 - t108 * t110) * t100, 0, 0, 0, 0, 0, -qJD(3) * t30 - qJDD(3) * t54 + (-t104 * t143 + t107 * t131) * t100, -qJD(3) * t29 - qJDD(3) * t55 + (-t104 * t131 - t107 * t143) * t100, 0, 0, 0, 0, 0, t23 * t54 + t30 * t67 + t31 * t64 - t6 * t88, t126 * t64 - t22 * t54 + t30 * t69 + t7 * t88, t126 * t23 + t22 * t31 - t6 * t69 - t67 * t7, t1 * t31 + t10 * t7 - t126 * t2 + t21 * t30 + t3 * t54 + t4 * t6 - g(3); 0, qJDD(2), t170 * t175 + t135, -t170 * t177 + t134, qJDD(2) * t97 + 0.2e1 * t104 * t144, 0.2e1 * t104 * t95 - 0.2e1 * t156 * t191, qJDD(3) * t104 + t107 * t109, qJDD(3) * t107 - t104 * t109, 0, t113 * t104 + t107 * t211, -t104 * t211 + t113 * t107, t107 * t69 * t159 + (-t106 * t22 - t161 * t69) * t104, (-t103 * t69 - t106 * t67) * t162 + (t185 - t106 * t23 + (t103 * t67 - t187) * qJD(4)) * t104, (t22 - t152) * t107 + (qJD(3) * t69 + t127) * t104, (t164 * t88 + t23) * t107 + (-t128 - t181) * t104, -t107 * t64 - t163 * t88, t66 * t64 + t216 * t88 + (t180 * t75 + t119) * t103 + (pkin(7) * t181 - t19 + (-pkin(7) * t64 + qJD(3) * t36 + qJD(4) * t44 + t12) * t103 + t209 * t106) * t107 + (pkin(7) * t23 + qJD(3) * t15 + t13 * t103 - t149 * t67 + t160 * t36) * t104, -t192 * t64 + t215 * t88 + t119 * t106 + ((pkin(7) * t69 + t106 * t36) * qJD(3) - t209 * t103 + t153) * t107 + (-t69 * t149 - t36 * t161 - t16 * qJD(3) + t13 * t106 + (-t22 - t152) * pkin(7)) * t104, t22 * t24 - t23 * t33 - t204 * t69 - t203 * t67 + (-t10 * t103 - t106 * t4) * t162 + (-t1 * t106 - t103 * t2 + (-t10 * t106 + t103 * t4) * qJD(4) - t118) * t104, t2 * t33 + t1 * t24 - g(1) * (t151 * t52 + t217 * t51) - g(2) * (t151 * t50 + t217 * t49) + t204 * t4 + t21 * t151 * t162 + (pkin(4) * t160 * t21 + t135 * t198 + t151 * t3) * t104 + t203 * t10 + (-g(3) * t151 * t105 + (-t21 * t168 - g(3) * (t104 * t198 - t217)) * t108) * t100; 0, 0, 0, 0, -t104 * t110 * t107, t191 * t110, t154, t95, qJDD(3), qJD(3) * t43 - t166 * t74 + t112, -t141 + (-qJD(2) * t74 - t46) * t107 + t121, -t187 * t88 - t185, (-t22 + t200) * t106 + (-t23 + t199) * t103, (-t104 * t69 + t172 * t88) * qJD(2) + t128, (t104 * t67 - t107 * t189) * qJD(2) + t127, t88 * t166, -t15 * t166 - pkin(3) * t23 - t43 * t67 + t58 * t88 + (-t214 * t88 + t117) * t103 + t114 * t106, pkin(3) * t22 - t103 * t114 + t106 * t117 + t16 * t166 - t197 * t88 - t43 * t69, -t22 * t76 - t23 * t77 - t195 * t69 - t196 * t67 + (t4 * t88 + t2) * t106 + (t10 * t88 - t1) * t103 - t121, t2 * t77 - t1 * t76 - t3 * t92 - g(1) * (t198 * t28 - t27 * t92) - g(2) * (t198 * t26 - t25 * t92) - g(3) * (t198 * t55 - t54 * t92) + t195 * t4 + (-pkin(4) * t189 - t43) * t21 + t196 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t67, -t63 + t208, -t22 - t200, -t199 - t23, t64, -t16 * t88 - t36 * t69 + t115 + t213, -t15 * t88 + t36 * t67 - g(1) * (-t103 * t51 - t106 * t28) - g(2) * (-t103 * t49 - t106 * t26) - g(3) * t126 + t123, pkin(4) * t22 - t207 * t67, t207 * t10 + (-t21 * t69 + t1 + t213) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63 - t208, t10 * t67 + t4 * t69 - t112 - t178 + t212;];
tau_reg = t5;
