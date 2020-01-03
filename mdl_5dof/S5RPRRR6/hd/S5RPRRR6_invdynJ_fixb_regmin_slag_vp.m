% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:43
% EndTime: 2019-12-31 19:01:48
% DurationCPUTime: 1.73s
% Computational Cost: add. (2174->237), mult. (4594->337), div. (0->0), fcn. (3286->14), ass. (0->155)
t109 = qJD(3) + qJD(4);
t117 = sin(qJ(4));
t121 = cos(qJ(3));
t199 = cos(qJ(4));
t166 = qJD(1) * t199;
t118 = sin(qJ(3));
t179 = qJD(1) * t118;
t208 = -t117 * t179 + t121 * t166;
t210 = t208 * t109;
t113 = qJ(3) + qJ(4);
t106 = sin(t113);
t110 = qJ(1) + pkin(9);
t102 = sin(t110);
t103 = cos(t110);
t206 = g(1) * t103 + g(2) * t102;
t209 = t206 * t106;
t114 = sin(pkin(9));
t96 = t114 * pkin(1) + pkin(6);
t200 = pkin(7) + t96;
t100 = t117 * pkin(3) + pkin(8);
t182 = t117 * t121;
t63 = -qJD(1) * t182 - t118 * t166;
t43 = -t63 * pkin(4) - pkin(8) * t208;
t60 = qJD(5) - t208;
t205 = (pkin(3) * t179 + qJD(5) * t100 + t43) * t60;
t204 = (pkin(8) * qJD(5) + t43) * t60;
t162 = t200 * qJD(1);
t51 = t121 * qJD(2) - t162 * t118;
t52 = t118 * qJD(2) + t162 * t121;
t170 = t199 * t52;
t187 = qJD(3) * pkin(3);
t47 = t51 + t187;
t26 = t117 * t47 + t170;
t104 = t121 * qJDD(2);
t78 = t96 * qJDD(1);
t155 = pkin(7) * qJDD(1) + t78;
t28 = qJDD(3) * pkin(3) - t52 * qJD(3) - t155 * t118 + t104;
t33 = t51 * qJD(3) + t118 * qJDD(2) + t155 * t121;
t203 = t26 * qJD(4) + t117 * t33 - t199 * t28;
t108 = qJDD(3) + qJDD(4);
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t156 = qJDD(1) * t199;
t172 = t121 * qJDD(1);
t35 = t117 * t172 + t118 * t156 + t210;
t50 = t116 * t109 - t120 * t63;
t14 = t50 * qJD(5) - t120 * t108 + t116 * t35;
t173 = t118 * qJDD(1);
t144 = t117 * t173 - t121 * t156;
t70 = t199 * t118 + t182;
t45 = t109 * t70;
t36 = t45 * qJD(1) + t144;
t32 = qJDD(5) + t36;
t136 = -t117 * t118 + t199 * t121;
t44 = t109 * t136;
t145 = t32 * t70 + t44 * t60;
t177 = qJD(5) * t116;
t169 = t70 * t177;
t202 = -t120 * t145 + t60 * t169;
t65 = t200 * t118;
t66 = t200 * t121;
t138 = -t117 * t66 - t199 * t65;
t165 = qJD(3) * t200;
t58 = t118 * t165;
t59 = t121 * t165;
t15 = t138 * qJD(4) - t117 * t59 - t199 * t58;
t164 = t199 * qJD(4);
t178 = qJD(4) * t117;
t129 = t117 * t28 + t47 * t164 - t52 * t178 + t199 * t33;
t115 = cos(pkin(9));
t97 = -t115 * pkin(1) - pkin(2);
t75 = -t121 * pkin(3) + t97;
t64 = t75 * qJD(1);
t37 = -pkin(4) * t208 + t63 * pkin(8) + t64;
t159 = t108 * pkin(8) + qJD(5) * t37 + t129;
t186 = t117 * t52;
t25 = t199 * t47 - t186;
t22 = -t109 * pkin(4) - t25;
t4 = -t108 * pkin(4) + t203;
t40 = -pkin(4) * t136 - t70 * pkin(8) + t75;
t42 = -t117 * t65 + t199 * t66;
t201 = -(qJD(5) * t40 + t15) * t60 + t159 * t136 + t22 * t44 - t42 * t32 + t4 * t70;
t98 = g(3) * t106;
t107 = cos(t113);
t196 = g(3) * t107;
t195 = t22 * t208;
t194 = t22 * t70;
t193 = t40 * t32;
t48 = -t120 * t109 - t116 * t63;
t192 = t48 * t60;
t191 = t50 * t60;
t190 = t60 * t63;
t189 = t63 * t208;
t176 = qJD(5) * t120;
t13 = t116 * t108 + t109 * t176 + t120 * t35 + t63 * t177;
t188 = -t13 * t136 + t50 * t45;
t185 = t13 * t116;
t81 = qJD(1) * t97;
t184 = t107 * t116;
t183 = t107 * t120;
t181 = qJDD(2) - g(3);
t111 = t118 ^ 2;
t180 = -t121 ^ 2 + t111;
t174 = qJD(1) * qJD(3);
t171 = t118 * t187;
t168 = -t4 - t196;
t163 = t118 * t174;
t23 = t109 * pkin(8) + t26;
t53 = pkin(3) * t163 + t75 * qJDD(1);
t9 = t36 * pkin(4) - t35 * pkin(8) + t53;
t160 = qJD(5) * t23 - t9;
t153 = t120 * t60;
t29 = t117 * t51 + t170;
t150 = pkin(3) * t178 - t29;
t148 = g(1) * t102 - g(2) * t103;
t119 = sin(qJ(1));
t122 = cos(qJ(1));
t147 = g(1) * t119 - g(2) * t122;
t146 = t136 * t14 - t45 * t48;
t11 = t116 * t37 + t120 * t23;
t143 = g(3) * t184 - t11 * t63 + t4 * t116 + t22 * t176;
t141 = t70 * t108 + t44 * t109;
t10 = -t116 * t23 + t120 * t37;
t140 = t10 * t63 + t120 * t209 + t22 * t177;
t139 = -t159 + t98;
t135 = -pkin(8) * t32 + t25 * t60 - t195;
t133 = -t81 * qJD(1) + t206 - t78;
t132 = 0.2e1 * t81 * qJD(3) - qJDD(3) * t96;
t123 = qJD(3) ^ 2;
t131 = -0.2e1 * qJDD(1) * t97 - t123 * t96 + t148;
t130 = -t70 * t60 * t176 - t145 * t116;
t30 = t199 * t51 - t186;
t128 = -t100 * t32 - t195 + (-pkin(3) * t164 + t30) * t60;
t127 = t206 * t107 - t208 * t64 - t129 + t98;
t126 = t64 * t63 - t196 - t203 + t209;
t124 = qJD(1) ^ 2;
t101 = -t199 * pkin(3) - pkin(4);
t77 = qJDD(3) * t121 - t123 * t118;
t76 = qJDD(3) * t118 + t123 * t121;
t57 = t102 * t116 + t103 * t183;
t56 = t102 * t120 - t103 * t184;
t55 = -t102 * t183 + t103 * t116;
t54 = t102 * t184 + t103 * t120;
t38 = -t208 ^ 2 + t63 ^ 2;
t34 = t108 * t136 - t45 * t109;
t19 = t45 * pkin(4) - t44 * pkin(8) + t171;
t18 = -t144 + (-qJD(1) * t70 - t63) * t109;
t17 = t35 - t210;
t16 = t42 * qJD(4) - t117 * t58 + t199 * t59;
t8 = t120 * t9;
t7 = t116 * t32 + t60 * t153 + t50 * t63;
t6 = -t60 ^ 2 * t116 + t120 * t32 - t48 * t63;
t5 = t50 * t153 + t185;
t1 = (t13 - t192) * t120 + (-t14 - t191) * t116;
t2 = [qJDD(1), t147, g(1) * t122 + g(2) * t119, (t147 + (t114 ^ 2 + t115 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t111 * qJDD(1) + 0.2e1 * t121 * t163, 0.2e1 * t118 * t172 - 0.2e1 * t180 * t174, t76, t77, 0, t132 * t118 + t131 * t121, -t131 * t118 + t132 * t121, t35 * t70 - t63 * t44, t136 * t35 + t208 * t44 - t70 * t36 + t63 * t45, t141, t34, 0, t148 * t107 + t108 * t138 - t16 * t109 - t136 * t53 - t171 * t208 + t75 * t36 + t64 * t45, -t148 * t106 - t42 * t108 - t15 * t109 - t63 * t171 + t75 * t35 + t64 * t44 + t53 * t70, -t50 * t169 + (t13 * t70 + t44 * t50) * t120, (-t116 * t50 - t120 * t48) * t44 + (-t185 - t120 * t14 + (t116 * t48 - t120 * t50) * qJD(5)) * t70, t188 - t202, t130 + t146, -t136 * t32 + t60 * t45, -g(1) * t55 - g(2) * t57 + t10 * t45 - t138 * t14 + t16 * t48 - t8 * t136 + (t19 * t60 + t193 + (t136 * t23 - t42 * t60 + t194) * qJD(5)) * t120 + t201 * t116, -g(1) * t54 - g(2) * t56 - t11 * t45 - t138 * t13 + t16 * t50 + (-(-qJD(5) * t42 + t19) * t60 - t193 - t160 * t136 - qJD(5) * t194) * t116 + t201 * t120; 0, 0, 0, t181, 0, 0, 0, 0, 0, t77, -t76, 0, 0, 0, 0, 0, t34, -t141, 0, 0, 0, 0, 0, t130 - t146, t188 + t202; 0, 0, 0, 0, -t118 * t124 * t121, t180 * t124, t173, t172, qJDD(3), -g(3) * t121 + t133 * t118 + t104, -t181 * t118 + t133 * t121, t189, t38, t17, t18, t108, t29 * t109 + (t199 * t108 - t109 * t178 + t179 * t208) * pkin(3) + t126, t30 * t109 + (-t108 * t117 - t109 * t164 + t63 * t179) * pkin(3) + t127, t5, t1, t7, t6, t190, t101 * t14 + t150 * t48 + (t168 - t205) * t120 + t128 * t116 + t140, t101 * t13 + t150 * t50 + t128 * t120 + (-t209 + t205) * t116 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t38, t17, t18, t108, t26 * t109 + t126, t25 * t109 + t127, t5, t1, t7, t6, t190, -pkin(4) * t14 - t26 * t48 + t135 * t116 + (t168 - t204) * t120 + t140, -pkin(4) * t13 - t26 * t50 + t135 * t120 + (-t209 + t204) * t116 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t48, -t48 ^ 2 + t50 ^ 2, t13 + t192, -t14 + t191, t32, -g(1) * t56 + g(2) * t54 + t11 * t60 + t116 * t139 - t23 * t176 - t22 * t50 + t8, g(1) * t57 - g(2) * t55 + t10 * t60 + t116 * t160 + t120 * t139 + t22 * t48;];
tau_reg = t2;
