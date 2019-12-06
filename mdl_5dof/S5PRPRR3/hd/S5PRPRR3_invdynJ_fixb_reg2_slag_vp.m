% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRR3
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:43
% EndTime: 2019-12-05 15:47:49
% DurationCPUTime: 1.81s
% Computational Cost: add. (2170->261), mult. (4540->348), div. (0->0), fcn. (3393->14), ass. (0->158)
t124 = qJD(4) ^ 2;
t114 = sin(pkin(8));
t116 = cos(pkin(8));
t150 = g(1) * t116 + g(2) * t114;
t109 = qJ(2) + pkin(9);
t97 = sin(t109);
t98 = cos(t109);
t136 = -g(3) * t98 + t150 * t97;
t113 = sin(pkin(9));
t115 = cos(pkin(9));
t119 = sin(qJ(2));
t122 = cos(qJ(2));
t68 = t113 * t122 + t115 * t119;
t58 = t68 * qJD(1);
t130 = qJD(2) * t58 + t136;
t172 = qJD(1) * qJD(2);
t213 = qJDD(1) * t119 + t122 * t172;
t101 = t122 * qJDD(1);
t191 = qJDD(2) * pkin(2);
t64 = -t119 * t172 + t101 + t191;
t165 = t213 * t113 - t115 * t64;
t30 = -qJDD(2) * pkin(3) + t165;
t92 = pkin(2) * t113 + pkin(6);
t208 = pkin(2) * t115;
t93 = -pkin(3) - t208;
t216 = -qJDD(2) * t93 - t124 * t92 + t130 - t30;
t117 = sin(qJ(5));
t120 = cos(qJ(5));
t121 = cos(qJ(4));
t167 = t121 * qJDD(2);
t118 = sin(qJ(4));
t168 = t118 * qJDD(2);
t108 = qJD(4) + qJD(5);
t70 = t117 * t121 + t118 * t120;
t215 = t108 * t70;
t22 = qJD(2) * t215 + t117 * t168 - t120 * t167;
t37 = t113 * t64 + t213 * t115;
t31 = qJDD(2) * pkin(6) + t37;
t214 = qJD(4) * qJD(3) + t31;
t185 = t120 * t121;
t186 = t117 * t118;
t69 = -t185 + t186;
t33 = t69 * t68;
t156 = -g(1) * t114 + g(2) * t116;
t211 = g(3) * t97;
t137 = t150 * t98 + t211;
t175 = qJD(5) * t117;
t102 = t121 * qJD(3);
t180 = qJD(1) * t119;
t179 = qJD(1) * t122;
t81 = qJD(2) * pkin(2) + t179;
t52 = t113 * t81 + t115 * t180;
t46 = qJD(2) * pkin(6) + t52;
t157 = pkin(7) * qJD(2) + t46;
t34 = -t157 * t118 + t102;
t29 = qJD(4) * pkin(4) + t34;
t177 = qJD(3) * t118;
t35 = t157 * t121 + t177;
t100 = t121 * qJDD(3);
t6 = qJDD(4) * pkin(4) + t100 + (-pkin(7) * qJDD(2) - t31) * t118 - t35 * qJD(4);
t166 = t118 * qJDD(3) + t214 * t121;
t176 = qJD(4) * t118;
t12 = -t46 * t176 + t166;
t171 = qJD(2) * qJD(4);
t141 = t118 * t171 - t167;
t7 = -pkin(7) * t141 + t12;
t1 = (qJD(5) * t29 + t7) * t120 + t117 * t6 - t35 * t175;
t209 = pkin(7) + t92;
t207 = pkin(2) * t119;
t162 = qJD(2) * t185;
t178 = qJD(2) * t118;
t163 = t117 * t178;
t61 = -t162 + t163;
t63 = t70 * qJD(2);
t203 = t63 * t61;
t65 = t209 * t118;
t66 = t209 * t121;
t38 = -t117 * t66 - t120 * t65;
t160 = qJD(4) * t209;
t55 = t118 * t160;
t56 = t121 * t160;
t84 = t113 * t180;
t60 = t115 * t179 - t84;
t202 = qJD(5) * t38 - t117 * t56 - t120 * t55 + t69 * t60;
t39 = -t117 * t65 + t120 * t66;
t201 = -qJD(5) * t39 + t117 * t55 - t120 * t56 + t70 * t60;
t146 = t108 * t186;
t174 = qJD(5) * t120;
t42 = -qJD(4) * t185 - t121 * t174 + t146;
t200 = -t70 * t22 + t42 * t61;
t199 = t117 * t35;
t198 = t118 * t46;
t197 = t120 * t35;
t67 = t113 * t119 - t115 * t122;
t192 = qJD(2) * t67;
t194 = qJD(2) * t192;
t193 = qJD(2) * t60;
t112 = qJ(4) + qJ(5);
t103 = sin(t112);
t190 = t103 * t114;
t189 = t103 * t116;
t104 = cos(t112);
t188 = t104 * t114;
t187 = t104 * t116;
t41 = t121 * t46 + t177;
t184 = t41 * qJD(4);
t183 = qJDD(1) - g(3);
t110 = t118 ^ 2;
t111 = t121 ^ 2;
t182 = t110 - t111;
t181 = t110 + t111;
t125 = qJD(2) ^ 2;
t164 = t118 * t125 * t121;
t96 = pkin(4) * t121 + pkin(3);
t158 = t121 * t171;
t51 = t115 * t81 - t84;
t154 = -qJD(5) * t162 - t117 * t167 + (-t158 - t168) * t120;
t153 = qJDD(2) * t181;
t152 = t118 * t158;
t151 = pkin(4) * t176 - t58;
t21 = qJD(2) * t146 + t154;
t148 = -t21 * t69 + t215 * t63;
t107 = qJDD(4) + qJDD(5);
t145 = t107 * t70 - t108 * t42;
t9 = t117 * t29 + t197;
t40 = t102 - t198;
t144 = t118 * t40 - t121 * t41;
t57 = t68 * qJD(2);
t143 = -qJD(2) * t57 - qJDD(2) * t67;
t139 = t124 * t68 - t143;
t138 = 0.2e1 * t192 * qJD(4) - qJDD(4) * t68;
t133 = -g(3) * t122 + t150 * t119;
t45 = -qJD(2) * pkin(3) - t51;
t132 = -qJDD(4) * t92 + (qJD(2) * t93 + t45 + t60) * qJD(4);
t2 = -qJD(5) * t9 - t117 * t7 + t120 * t6;
t131 = -t45 * qJD(2) + t137;
t13 = -t118 * t31 + t100 - t184;
t129 = -t13 * t118 + t12 * t121 + (-t118 * t41 - t121 * t40) * qJD(4);
t44 = -t96 * qJD(2) - t51;
t127 = -g(1) * (-t98 * t187 - t190) - g(2) * (-t98 * t188 + t189) + t44 * t61 + t104 * t211 - t1;
t126 = -g(1) * (-t98 * t189 + t188) - g(2) * (-t98 * t190 - t187) - t44 * t63 + t2 + t103 * t211;
t123 = -pkin(7) - pkin(6);
t106 = t122 * pkin(2);
t76 = qJDD(4) * t121 - t118 * t124;
t75 = qJDD(4) * t118 + t121 * t124;
t74 = -t96 - t208;
t32 = t70 * t68;
t23 = -t61 ^ 2 + t63 ^ 2;
t20 = -t107 * t69 - t108 * t215;
t19 = pkin(4) * t141 + t30;
t17 = t108 * t63 - t22;
t16 = -t154 + (-t163 + t61) * t108;
t11 = t120 * t34 - t199;
t10 = -t117 * t34 - t197;
t8 = t120 * t29 - t199;
t4 = t108 * t33 + t192 * t70;
t3 = t192 * t69 - t215 * t68;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t183, 0, 0, 0, 0, 0, 0, qJDD(2) * t122 - t119 * t125, -qJDD(2) * t119 - t122 * t125, 0, -g(3) + (t119 ^ 2 + t122 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t143, -qJDD(2) * t68 + t194, 0, t165 * t67 - t192 * t52 + t37 * t68 - t51 * t57 - g(3), 0, 0, 0, 0, 0, 0, t118 * t138 - t121 * t139, t118 * t139 + t121 * t138, t68 * t153 - t181 * t194, t129 * t68 + t144 * t192 + t30 * t67 + t45 * t57 - g(3), 0, 0, 0, 0, 0, 0, -t107 * t32 + t108 * t4 + t22 * t67 + t57 * t61, t107 * t33 - t108 * t3 - t21 * t67 + t57 * t63, -t21 * t32 + t22 * t33 - t3 * t61 - t4 * t63, -t1 * t33 + t19 * t67 - t2 * t32 + t3 * t9 + t4 * t8 + t44 * t57 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t101 + t133, -t183 * t119 + t150 * t122, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t115 * t191 + t130 - t165, -t113 * t191 + t137 + t193 - t37, 0, t51 * t58 - t52 * t60 + (t113 * t37 - t115 * t165 + t133) * pkin(2), qJDD(2) * t110 + 0.2e1 * t152, 0.2e1 * t118 * t167 - 0.2e1 * t182 * t171, t75, qJDD(2) * t111 - 0.2e1 * t152, t76, 0, t132 * t118 + t216 * t121, -t216 * t118 + t132 * t121, t92 * t153 - t181 * t193 + t129 - t137, t30 * t93 - t45 * t58 - g(3) * (pkin(3) * t98 + pkin(6) * t97 + t106) + t144 * t60 + t129 * t92 + t150 * (pkin(3) * t97 - pkin(6) * t98 + t207), -t21 * t70 - t42 * t63, -t148 + t200, t145, t215 * t61 + t22 * t69, t20, 0, t136 * t104 + t107 * t38 + t201 * t108 + t151 * t61 + t19 * t69 + t215 * t44 + t22 * t74, -t103 * t136 - t107 * t39 - t202 * t108 + t151 * t63 + t19 * t70 - t21 * t74 - t42 * t44, -t1 * t69 - t2 * t70 - t201 * t63 - t202 * t61 + t21 * t38 - t215 * t9 - t22 * t39 + t42 * t8 - t137, t1 * t39 + t2 * t38 + t19 * t74 - g(3) * (-t123 * t97 + t96 * t98 + t106) + t202 * t9 + t201 * t8 + t151 * t44 + t150 * (t123 * t98 + t96 * t97 + t207); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + t156, 0, 0, 0, 0, 0, 0, t76, -t75, 0, -qJD(4) * t144 + t118 * t12 + t121 * t13 + t156, 0, 0, 0, 0, 0, 0, t20, -t145, t148 + t200, t1 * t70 - t2 * t69 - t215 * t8 - t42 * t9 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t182 * t125, t168, t164, t167, qJDD(4), t184 + t100 + (-qJD(4) * t46 + t156) * t121 + (t131 - t214) * t118, -t156 * t118 + t131 * t121 + (t40 + t198) * qJD(4) - t166, 0, 0, t203, t23, t16, -t203, t17, t107, -t10 * t108 + (t107 * t120 - t108 * t175 - t178 * t61) * pkin(4) + t126, t108 * t11 + (-t107 * t117 - t108 * t174 - t178 * t63) * pkin(4) + t127, (t10 + t9) * t63 + (t11 - t8) * t61 + (-t117 * t22 + t120 * t21 + (t117 * t63 - t120 * t61) * qJD(5)) * pkin(4), -t8 * t10 - t9 * t11 + (t1 * t117 + t2 * t120 + t156 * t121 + (-t44 * qJD(2) + t137) * t118 + (-t117 * t8 + t120 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, t23, t16, -t203, t17, t107, t108 * t9 + t126, t108 * t8 + t127, 0, 0;];
tau_reg = t5;
