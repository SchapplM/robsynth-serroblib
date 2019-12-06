% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:29
% EndTime: 2019-12-05 18:20:33
% DurationCPUTime: 1.47s
% Computational Cost: add. (2547->264), mult. (4039->354), div. (0->0), fcn. (2482->14), ass. (0->171)
t108 = qJ(1) + qJ(2);
t98 = pkin(8) + t108;
t88 = sin(t98);
t89 = cos(t98);
t210 = g(2) * t89 + g(3) * t88;
t117 = cos(qJ(2));
t192 = pkin(1) * qJD(2);
t157 = qJD(1) * t192;
t114 = sin(qJ(2));
t164 = qJDD(1) * t114;
t211 = pkin(1) * t164 + t117 * t157;
t100 = cos(t108);
t99 = sin(t108);
t209 = g(2) * t100 + g(3) * t99;
t105 = qJD(1) + qJD(2);
t109 = sin(pkin(9));
t103 = t109 ^ 2;
t111 = cos(pkin(9));
t104 = t111 ^ 2;
t169 = t103 + t104;
t208 = t105 * t169;
t112 = cos(pkin(8));
t193 = pkin(1) * qJD(1);
t158 = t117 * t193;
t110 = sin(pkin(8));
t159 = t114 * t193;
t76 = t110 * t159;
t55 = t112 * t158 - t76;
t170 = qJD(4) - t55;
t207 = pkin(2) * t99;
t206 = pkin(1) * t114;
t115 = sin(qJ(1));
t205 = pkin(1) * t115;
t204 = pkin(1) * t117;
t118 = cos(qJ(1));
t203 = pkin(1) * t118;
t202 = pkin(2) * t100;
t201 = pkin(2) * t110;
t200 = pkin(2) * t112;
t113 = sin(qJ(5));
t116 = cos(qJ(5));
t167 = qJD(4) * t111;
t173 = t111 * t116;
t174 = t111 * t113;
t131 = -pkin(4) * t111 - pkin(7) * t109 - pkin(3);
t59 = t131 - t200;
t87 = qJ(4) + t201;
t34 = t116 * t59 - t174 * t87;
t182 = qJD(5) * t34;
t172 = t112 * t114;
t126 = pkin(1) * (t110 * t117 + t172);
t53 = qJD(1) * t126;
t199 = -t113 * t53 + t116 * t167 - t173 * t55 + t182;
t35 = t113 * t59 + t173 * t87;
t181 = qJD(5) * t35;
t198 = -t113 * t167 - t116 * t53 + t174 * t55 - t181;
t197 = t210 * t109;
t185 = t111 * t89;
t186 = t111 * t88;
t196 = g(2) * t185 + g(3) * t186;
t194 = g(2) * t88 - g(3) * t89;
t94 = pkin(2) + t204;
t58 = pkin(1) * t172 + t110 * t94;
t62 = pkin(2) * t105 + t158;
t41 = t112 * t62 - t76;
t139 = qJD(4) - t41;
t23 = t105 * t131 + t139;
t42 = t110 * t62 + t112 * t159;
t37 = qJ(4) * t105 + t42;
t32 = qJD(3) * t109 + t111 * t37;
t8 = t113 * t23 + t116 * t32;
t191 = qJD(5) * t8;
t190 = t105 * t53;
t54 = qJD(2) * t126;
t189 = t105 * t54;
t102 = qJDD(1) + qJDD(2);
t95 = qJDD(1) * t204;
t50 = pkin(2) * t102 - t114 * t157 + t95;
t30 = t110 * t50 + t211 * t112;
t20 = qJ(4) * t102 + qJD(4) * t105 + t30;
t14 = -t111 * qJDD(3) + t109 * t20;
t188 = t109 * t14;
t31 = -t111 * qJD(3) + t109 * t37;
t187 = t109 * t31;
t175 = t111 * t102;
t70 = -qJDD(5) + t175;
t184 = t70 * t111;
t85 = t110 * t206;
t57 = t112 * t94 - t85;
t38 = t131 - t57;
t51 = qJ(4) + t58;
t18 = t116 * t38 - t174 * t51;
t183 = qJD(5) * t18;
t180 = t102 * t113;
t179 = t102 * t116;
t101 = t105 ^ 2;
t178 = t103 * t101;
t177 = t105 * t111;
t176 = t109 * t102;
t171 = t113 * t116;
t106 = t113 ^ 2;
t107 = t116 ^ 2;
t168 = t106 - t107;
t166 = qJD(5) * t113;
t165 = qJD(5) * t116;
t153 = t109 * t166;
t7 = -t113 * t32 + t116 * t23;
t163 = t7 * t153 + t197;
t161 = t211 * t110 - t112 * t50;
t145 = qJDD(4) + t161;
t22 = -pkin(3) * t102 + t145;
t162 = t22 * t109 - t197;
t160 = t95 + t209;
t73 = -qJD(5) + t177;
t155 = t73 * t166;
t154 = -g(2) * t99 + g(3) * t100;
t152 = -t111 * t14 - g(1);
t151 = t70 - t175;
t150 = t70 + t175;
t149 = t102 * t169;
t148 = t105 * (-qJD(5) - t73);
t147 = qJD(1) * (-qJD(2) + t105);
t146 = qJD(2) * (-qJD(1) - t105);
t143 = t171 * t178;
t142 = t105 * t113 * t165;
t141 = t161 - t210;
t140 = -t30 - t194;
t138 = g(2) * t118 + g(3) * t115;
t137 = t113 * t7 - t116 * t8;
t56 = t112 * t117 * t192 - qJD(2) * t85;
t49 = qJD(4) + t56;
t136 = t14 * t51 + t31 * t49;
t52 = -pkin(3) - t57;
t135 = t102 * t52 + t189;
t90 = -pkin(3) - t200;
t134 = t102 * t90 - t190;
t15 = qJDD(3) * t109 + t111 * t20;
t133 = t15 * t111 + t188 + t194;
t132 = -pkin(3) * t88 + t89 * qJ(4) - t207;
t130 = t210 * pkin(7);
t129 = g(1) * t109 - qJD(5) * t23 - t15;
t19 = t113 * t38 + t173 * t51;
t128 = -qJD(5) * t32 - t105 * t187;
t127 = t14 * t87 + t170 * t31;
t125 = -pkin(3) * t89 - qJ(4) * t88 - t202;
t13 = t102 * t131 + t145;
t11 = t116 * t13;
t3 = -t113 * t15 + t11 - t191;
t44 = -t89 * t113 + t173 * t88;
t46 = -t113 * t88 - t173 * t89;
t124 = -g(2) * t46 + g(3) * t44 - t111 * t3 + t113 * t188 + t165 * t187;
t123 = -t73 ^ 2 - t178;
t122 = -pkin(4) * t186 + t132;
t121 = -pkin(4) * t185 + t125;
t2 = qJD(5) * t7 + t113 * t13 + t116 * t15;
t43 = t116 * t89 + t174 * t88;
t45 = -t116 * t88 + t174 * t89;
t120 = -g(2) * t45 - g(3) * t43 + t2 * t111 + t116 * t188 - t153 * t31;
t84 = t104 * t102;
t83 = t103 * t102;
t63 = 0.2e1 * t109 * t175;
t61 = t153 * t177;
t40 = (t102 * t107 - 0.2e1 * t142) * t103;
t39 = (t102 * t106 + 0.2e1 * t142) * t103;
t36 = -pkin(3) * t105 + t139;
t33 = 0.2e1 * (qJD(5) * t105 * t168 - t102 * t171) * t103;
t17 = (t150 * t113 + (t73 + t177) * t165) * t109;
t16 = t61 + (-t116 * t150 + t155) * t109;
t5 = -qJD(5) * t19 + t116 * t54 - t174 * t49;
t4 = t113 * t54 + t173 * t49 + t183;
t1 = [0, 0, 0, 0, 0, qJDD(1), t138, -g(2) * t115 + g(3) * t118, 0, 0, 0, 0, 0, 0, 0, t102, (t102 * t117 + t114 * t146) * pkin(1) + t160, ((-qJDD(1) - t102) * t114 + t117 * t146) * pkin(1) + t154, 0, (t138 + (t114 ^ 2 + t117 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t102, t102 * t57 - t141 - t189, -t102 * t58 - t105 * t56 + t140, 0, t30 * t58 + t42 * t56 - t161 * t57 - t41 * t54 - g(2) * (-t202 - t203) - g(3) * (-t205 - t207), t83, t63, 0, t84, 0, 0, (-t135 - t22) * t111 + t196, t109 * t135 + t162, t149 * t51 + t49 * t208 + t133, t22 * t52 + t36 * t54 - g(2) * (t125 - t203) - g(3) * (t132 - t205) + (t15 * t51 + t32 * t49) * t111 + t136 * t109, t40, t33, t16, t39, t17, t184, -t18 * t70 - t5 * t73 + (t51 * t180 + (t113 * t49 + t165 * t51) * t105) * t103 + t124, t19 * t70 + t4 * t73 + (t51 * t179 + (t116 * t49 - t166 * t51) * t105) * t103 + t120, ((-t102 * t19 - t2 + (-t4 + t183) * t105) * t113 + (-t102 * t18 - t105 * t5 - t3 + (-t105 * t19 - t8) * qJD(5)) * t116) * t109 + t163, t2 * t19 + t8 * t4 + t3 * t18 + t7 * t5 - g(2) * (t121 - t203) - g(3) * (t122 - t205) + (t130 + t136) * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t147 * t206 + t160, (t117 * t147 - t164) * pkin(1) + t154, 0, 0, 0, 0, 0, 0, 0, t102, t102 * t200 - t141 + t190, -t102 * t201 + t105 * t55 + t140, 0, t41 * t53 - t42 * t55 + (t110 * t30 - t112 * t161 + t209) * pkin(2), t83, t63, 0, t84, 0, 0, (-t134 - t22) * t111 + t196, t109 * t134 + t162, t87 * t149 + t170 * t208 + t133, t22 * t90 - t36 * t53 - g(2) * t125 - g(3) * t132 + (t15 * t87 + t170 * t32) * t111 + t127 * t109, t40, t33, t16, t39, t17, t184, -t34 * t70 - t198 * t73 + (t87 * t180 + (t113 * t170 + t165 * t87) * t105) * t103 + t124, t35 * t70 + t199 * t73 + (t87 * t179 + (t116 * t170 - t166 * t87) * t105) * t103 + t120, ((-t102 * t35 - t2) * t113 + (-t102 * t34 - t191 - t3) * t116 + ((-t181 - t198) * t116 + (t182 - t199) * t113) * t105) * t109 + t163, t2 * t35 + t3 * t34 - g(2) * t121 - g(3) * t122 + t199 * t8 + t198 * t7 + (t130 + t127) * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t15 + t152, 0, 0, 0, 0, 0, 0, (t151 * t113 + (t73 - t177) * t165) * t109, t61 + (t116 * t151 - t155) * t109, 0, (-t113 * t3 + t116 * t2 + (-t113 * t8 - t116 * t7) * qJD(5)) * t109 + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, t176, -t169 * t101, (-t111 * t32 - t187) * t105 + t22 - t210, 0, 0, 0, 0, 0, 0, t113 * t123 - t116 * t70, t113 * t70 + t116 * t123, (-t106 - t107) * t176, t113 * t2 + t116 * t3 - t137 * qJD(5) + (t111 * t137 - t187) * t105 - t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, -t168 * t178, (t113 * t148 + t179) * t109, -t143, (t116 * t148 - t180) * t109, -t70, -g(2) * t43 + g(3) * t45 + t113 * t129 + t116 * t128 - t73 * t8 + t11, -g(2) * t44 - g(3) * t46 - t7 * t73 + t129 * t116 + (-t128 - t13) * t113, 0, 0;];
tau_reg = t1;
