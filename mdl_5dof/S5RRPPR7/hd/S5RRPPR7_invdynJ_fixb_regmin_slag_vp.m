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
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 19:36:25
% EndTime: 2019-12-31 19:36:30
% DurationCPUTime: 1.73s
% Computational Cost: add. (1867->290), mult. (4366->374), div. (0->0), fcn. (3077->10), ass. (0->154)
t107 = sin(pkin(8));
t108 = cos(pkin(8));
t111 = sin(qJ(2));
t114 = cos(qJ(2));
t78 = t107 * t114 + t108 * t111;
t70 = t78 * qJD(1);
t194 = qJD(5) + t70;
t110 = sin(qJ(5));
t113 = cos(qJ(5));
t161 = qJD(1) * t111;
t170 = t108 * t114;
t67 = -qJD(1) * t170 + t107 * t161;
t51 = qJD(2) * t110 - t113 * t67;
t195 = t194 * t51;
t145 = t110 * t194;
t157 = qJD(1) * qJD(2);
t151 = t114 * t157;
t152 = t111 * t157;
t45 = qJDD(1) * t78 - t107 * t152 + t108 * t151;
t40 = qJDD(5) + t45;
t133 = t113 * t40 - t145 * t194;
t66 = t70 ^ 2;
t192 = -t67 ^ 2 - t66;
t109 = -qJ(3) - pkin(6);
t84 = t109 * t114;
t81 = qJD(1) * t84;
t73 = t107 * t81;
t83 = t109 * t111;
t80 = qJD(1) * t83;
t48 = t108 * t80 + t73;
t165 = -qJD(4) + t48;
t191 = -qJD(5) + t194;
t112 = sin(qJ(1));
t115 = cos(qJ(1));
t190 = g(1) * t112 - g(2) * t115;
t142 = g(1) * t115 + g(2) * t112;
t185 = pkin(4) * t67;
t173 = t108 * t81;
t76 = qJD(2) * pkin(2) + t80;
t43 = t107 * t76 - t173;
t34 = -qJD(2) * qJ(4) - t43;
t21 = -t34 - t185;
t47 = t107 * t80 - t173;
t95 = -pkin(2) * t108 - pkin(3);
t90 = -pkin(7) + t95;
t189 = t90 * t40 + (t21 - t47 + t185) * t194;
t187 = pkin(3) + pkin(7);
t155 = t114 * qJDD(1);
t156 = t111 * qJDD(1);
t136 = -t107 * t156 + t108 * t155;
t69 = t78 * qJD(2);
t44 = qJD(1) * t69 - t136;
t186 = pkin(3) * t44;
t184 = pkin(4) * t70;
t183 = pkin(2) * t114;
t103 = qJ(2) + pkin(8);
t100 = cos(t103);
t93 = g(3) * t100;
t179 = g(3) * t114;
t77 = t107 * t111 - t170;
t178 = t21 * t77;
t96 = pkin(1) + t183;
t135 = -qJ(4) * t78 - t96;
t22 = t187 * t77 + t135;
t177 = t22 * t40;
t176 = t51 * t67;
t53 = qJD(2) * t113 + t110 * t67;
t175 = t53 * t67;
t146 = qJD(2) * t109;
t65 = -qJD(3) * t111 + t114 * t146;
t39 = qJDD(2) * pkin(2) + qJD(1) * t65 + qJDD(1) * t83;
t64 = qJD(3) * t114 + t111 * t146;
t46 = qJD(1) * t64 - qJDD(1) * t84;
t11 = t107 * t39 + t108 * t46;
t174 = qJ(4) * t45;
t158 = qJD(5) * t113;
t159 = qJD(5) * t110;
t13 = -qJD(2) * t159 + qJDD(2) * t113 + t110 * t44 + t158 * t67;
t172 = t113 * t13;
t171 = qJDD(2) * pkin(3);
t169 = t110 * t112;
t168 = t110 * t115;
t167 = t112 * t113;
t166 = t113 * t115;
t164 = t184 - t165;
t105 = t111 ^ 2;
t162 = -t114 ^ 2 + t105;
t160 = qJD(2) * t111;
t154 = qJDD(2) * qJ(4) + t11;
t98 = pkin(2) * t160;
t153 = t77 * t158;
t149 = pkin(2) * t161 + qJ(4) * t67;
t127 = pkin(2) * t152 - qJDD(1) * t96 + qJDD(3);
t121 = -qJD(4) * t70 + t127 - t174;
t1 = t187 * t44 + t121;
t42 = t108 * t76 + t73;
t143 = qJD(4) - t42;
t16 = -qJD(2) * t187 + t143 + t184;
t148 = qJD(5) * t16 + t1;
t82 = -qJD(1) * t96 + qJD(3);
t125 = -qJ(4) * t70 + t82;
t15 = t187 * t67 + t125;
t10 = -t107 * t46 + t108 * t39;
t138 = qJDD(4) - t10;
t5 = pkin(4) * t45 - qJDD(2) * t187 + t138;
t147 = -qJD(5) * t15 + t5;
t27 = t107 * t64 - t108 * t65;
t49 = -t107 * t84 - t108 * t83;
t144 = qJDD(2) * t110 - t113 * t44;
t139 = t194 * t69 + t40 * t77;
t99 = sin(t103);
t137 = pkin(3) * t100 + qJ(4) * t99;
t28 = t107 * t65 + t108 * t64;
t50 = t107 * t83 - t108 * t84;
t3 = t110 * t16 + t113 * t15;
t8 = -qJD(2) * qJD(4) - t154;
t72 = qJD(2) * t170 - t107 * t160;
t132 = -qJ(4) * t72 - qJD(4) * t78 + t98;
t131 = -0.2e1 * pkin(1) * t157 - pkin(6) * qJDD(2);
t29 = pkin(4) * t78 + t49;
t6 = -pkin(4) * t44 - t8;
t130 = t21 * t69 - t29 * t40 + t6 * t77;
t128 = -t113 * t194 ^ 2 - t110 * t40;
t126 = -g(3) * t99 - t100 * t142;
t116 = qJD(2) ^ 2;
t124 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t116 + t190;
t117 = qJD(1) ^ 2;
t123 = pkin(1) * t117 - pkin(6) * qJDD(1) + t142;
t122 = t127 - t190;
t120 = t27 * t70 - t28 * t67 - t44 * t50 + t45 * t49 - t142;
t25 = pkin(3) * t67 + t125;
t119 = -t142 * t99 + t25 * t70 + t138 + t93;
t118 = t6 + (-qJD(5) * t90 + t187 * t70 + t149) * t194 + t126;
t92 = pkin(2) * t107 + qJ(4);
t86 = t115 * t96;
t63 = -t169 * t99 + t166;
t62 = t167 * t99 + t168;
t61 = t168 * t99 + t167;
t60 = t166 * t99 - t169;
t58 = qJD(2) * t67;
t41 = pkin(3) * t77 + t135;
t31 = -qJD(2) * pkin(3) + t143;
t30 = -pkin(4) * t77 + t50;
t26 = pkin(3) * t70 + t149;
t20 = pkin(3) * t69 + t132;
t19 = -pkin(4) * t69 + t28;
t18 = pkin(4) * t72 + t27;
t14 = qJD(5) * t53 + t144;
t12 = t187 * t69 + t132;
t9 = t138 - t171;
t7 = t121 + t186;
t4 = t113 * t5;
t2 = -t110 * t15 + t113 * t16;
t17 = [qJDD(1), t190, t142, qJDD(1) * t105 + 0.2e1 * t111 * t151, 0.2e1 * t111 * t155 - 0.2e1 * t157 * t162, qJDD(2) * t111 + t114 * t116, qJDD(2) * t114 - t111 * t116, 0, t111 * t131 + t114 * t124, -t111 * t124 + t114 * t131, -t10 * t78 - t11 * t77 - t42 * t72 - t43 * t69 + t120, t11 * t50 + t43 * t28 - t10 * t49 - t42 * t27 - t127 * t96 + t82 * t98 - g(1) * (-t109 * t115 - t112 * t96) - g(2) * (-t109 * t112 + t86), t31 * t72 + t34 * t69 + t77 * t8 + t78 * t9 + t120, qJD(2) * t27 + qJDD(2) * t49 - t100 * t190 - t20 * t67 - t25 * t69 - t41 * t44 - t7 * t77, qJD(2) * t28 + qJDD(2) * t50 + t190 * t99 - t20 * t70 - t25 * t72 - t41 * t45 - t7 * t78, -g(2) * t86 + t25 * t20 + t31 * t27 - t34 * t28 + t7 * t41 + t9 * t49 - t8 * t50 + (g(1) * t109 - g(2) * t137) * t115 + (-g(1) * (-t137 - t96) + g(2) * t109) * t112, t53 * t153 + (t13 * t77 + t53 * t69) * t110, (-t110 * t51 + t113 * t53) * t69 + (-t110 * t14 + t172 + (-t110 * t53 - t113 * t51) * qJD(5)) * t77, t110 * t139 + t13 * t78 + t153 * t194 + t53 * t72, -t159 * t194 * t77 + t113 * t139 - t14 * t78 - t51 * t72, t194 * t72 + t40 * t78, -g(1) * t63 - g(2) * t61 + t30 * t14 + t19 * t51 + t2 * t72 + t4 * t78 + (-t1 * t78 - t12 * t194 - t177) * t110 + (t18 * t194 - t130) * t113 + ((-t110 * t29 - t113 * t22) * t194 - t3 * t78 + t110 * t178) * qJD(5), g(1) * t62 - g(2) * t60 + t30 * t13 + t19 * t53 - t3 * t72 + (-(qJD(5) * t29 + t12) * t194 - t177 - t148 * t78 + qJD(5) * t178) * t113 + (-(-qJD(5) * t22 + t18) * t194 - t147 * t78 + t130) * t110; 0, 0, 0, -t111 * t117 * t114, t162 * t117, t156, t155, qJDD(2), t111 * t123 - t179, g(3) * t111 + t114 * t123, (t43 - t47) * t70 + (-t42 + t48) * t67 + (-t107 * t44 - t108 * t45) * pkin(2), t42 * t47 - t43 * t48 + (-t179 + t10 * t108 + t107 * t11 + (-qJD(1) * t82 + t142) * t111) * pkin(2), -t44 * t92 + t45 * t95 + (-t34 - t47) * t70 + (t31 + t165) * t67, -qJD(2) * t47 + t26 * t67 + (-pkin(3) + t95) * qJDD(2) + t119, qJDD(2) * t92 - t25 * t67 + t26 * t70 + (0.2e1 * qJD(4) - t48) * qJD(2) + t126 + t154, -t8 * t92 + t9 * t95 - t25 * t26 - t31 * t47 - g(3) * (t137 + t183) + t165 * t34 + t142 * (pkin(2) * t111 + pkin(3) * t99 - qJ(4) * t100), -t145 * t53 + t172, (-t194 * t53 - t14) * t113 + (-t13 + t195) * t110, t133 + t175, t128 - t176, t194 * t67, t110 * t118 + t113 * t189 + t92 * t14 + t164 * t51 + t2 * t67, -t110 * t189 + t113 * t118 + t92 * t13 + t164 * t53 - t3 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t42 * t70 + t43 * t67 + t122, t192, -0.2e1 * qJD(2) * t70 + t136, -t45 + t58, t186 - t174 - t34 * t67 + (-qJD(4) - t31) * t70 + t122, 0, 0, 0, 0, 0, t128 + t176, t175 - t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 + t45, -t67 * t70 + qJDD(2), -t66 - t116, qJD(2) * t34 + t119 - t171, 0, 0, 0, 0, 0, -qJD(2) * t51 + t133, -qJD(2) * t53 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t51, -t51 ^ 2 + t53 ^ 2, t13 + t195, t191 * t53 - t144, t40, -g(1) * t60 - g(2) * t62 - t1 * t110 + t113 * t93 + t191 * t3 - t21 * t53 + t4, g(1) * t61 - g(2) * t63 + t2 * t194 + t21 * t51 - t148 * t113 + (-t147 - t93) * t110;];
tau_reg = t17;
