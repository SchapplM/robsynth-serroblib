% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:16
% EndTime: 2020-01-03 11:59:19
% DurationCPUTime: 1.52s
% Computational Cost: add. (1918->278), mult. (3139->325), div. (0->0), fcn. (1791->12), ass. (0->173)
t110 = qJDD(1) + qJDD(2);
t115 = sin(pkin(8));
t116 = cos(pkin(8));
t122 = cos(qJ(2));
t199 = pkin(1) * qJD(2);
t168 = qJD(1) * t199;
t119 = sin(qJ(2));
t175 = qJDD(1) * t119;
t212 = pkin(1) * t175 + t122 * t168;
t205 = t122 * pkin(1);
t97 = qJDD(1) * t205;
t44 = t110 * pkin(2) - t119 * t168 + t97;
t21 = t115 * t44 + t212 * t116;
t15 = t110 * pkin(7) + t21;
t213 = qJD(3) * qJD(4) + t15;
t118 = sin(qJ(4));
t121 = cos(qJ(4));
t102 = t121 * qJD(3);
t111 = qJD(1) + qJD(2);
t200 = pkin(1) * qJD(1);
t170 = t119 * t200;
t169 = t122 * t200;
t61 = t111 * pkin(2) + t169;
t35 = t115 * t61 + t116 * t170;
t30 = t111 * pkin(7) + t35;
t23 = -t118 * t30 + t102;
t177 = t118 * qJD(3);
t195 = t121 * t30;
t24 = t177 + t195;
t179 = qJD(4) * t118;
t5 = t118 * qJDD(3) + t213 * t121 - t30 * t179;
t99 = t121 * qJDD(3);
t6 = -t24 * qJD(4) - t118 * t15 + t99;
t127 = -t6 * t118 + t5 * t121 + (-t118 * t24 - t121 * t23) * qJD(4);
t114 = qJ(1) + qJ(2);
t100 = pkin(8) + t114;
t88 = cos(t100);
t81 = g(3) * t88;
t87 = sin(t100);
t82 = g(2) * t87;
t201 = t81 - t82;
t160 = qJ(5) * t111 + t30;
t139 = t160 * t121;
t211 = g(2) * t88;
t112 = t118 ^ 2;
t210 = pkin(4) * t112;
t209 = g(1) * t121;
t208 = t115 * pkin(2);
t207 = t116 * pkin(2);
t206 = t121 * pkin(4);
t16 = -t160 * t118 + t102;
t198 = qJD(4) * pkin(4);
t13 = t16 + t198;
t204 = t13 - t16;
t185 = t111 * t121;
t183 = t116 * t119;
t133 = pkin(1) * (t115 * t122 + t183);
t49 = qJD(1) * t133;
t77 = t115 * t170;
t51 = t116 * t169 - t77;
t203 = t51 * t179 + t49 * t185;
t196 = t118 * t87;
t202 = g(3) * t196 + t118 * t211;
t96 = pkin(2) + t205;
t54 = pkin(1) * t183 + t115 * t96;
t95 = pkin(3) + t206;
t62 = -t95 - t207;
t197 = t110 * t62;
t194 = t49 * t111;
t50 = qJD(2) * t133;
t193 = t50 * t111;
t184 = t115 * t119;
t52 = (t116 * t122 - t184) * t199;
t192 = t52 * t111;
t48 = pkin(7) + t54;
t191 = -qJ(5) - t48;
t89 = pkin(7) + t208;
t190 = -qJ(5) - t89;
t189 = qJ(5) * t110;
t188 = qJDD(4) * pkin(4);
t109 = t111 ^ 2;
t187 = t109 * t121;
t186 = t111 * t118;
t85 = t118 * t110;
t86 = t121 * t110;
t113 = t121 ^ 2;
t182 = t112 - t113;
t181 = t112 + t113;
t180 = qJD(4) * t111;
t178 = qJD(4) * t121;
t101 = t121 * qJD(5);
t166 = t111 * t179;
t3 = t111 * t101 + (-t166 + t86) * qJ(5) + t5;
t174 = t3 * t121 + t201;
t117 = -qJ(5) - pkin(7);
t103 = sin(t114);
t93 = pkin(2) * t103;
t173 = t88 * t117 + t87 * t95 + t93;
t172 = t212 * t115 - t116 * t44;
t104 = cos(t114);
t94 = pkin(2) * t104;
t171 = t88 * pkin(3) + t87 * pkin(7) + t94;
t165 = t111 * t178;
t164 = g(2) * t103 - g(3) * t104;
t34 = t116 * t61 - t77;
t22 = -t95 * t111 + qJD(5) - t34;
t74 = pkin(4) * t166;
t8 = -t95 * t110 + qJDD(5) + t172 + t74;
t163 = t8 * t118 + t22 * t178 + t202;
t53 = -pkin(1) * t184 + t116 * t96;
t162 = t181 * t51;
t14 = -t110 * pkin(3) + t172;
t29 = -t111 * pkin(3) - t34;
t161 = t14 * t118 + t29 * t178 + t202;
t159 = qJD(4) * t191;
t158 = qJD(4) * t190;
t157 = t110 * t181;
t155 = qJD(1) * (-qJD(2) + t111);
t154 = qJD(2) * (-qJD(1) - t111);
t153 = t87 * pkin(3) - t88 * pkin(7) + t93;
t47 = -pkin(3) - t53;
t151 = g(2) * t196 - t209 + t99;
t150 = t118 * t165;
t149 = -t87 * t117 + t88 * t95 + t94;
t148 = -t21 - t201;
t147 = g(3) * t87 + t211;
t146 = -g(2) * t104 - g(3) * t103;
t120 = sin(qJ(1));
t123 = cos(qJ(1));
t145 = -g(2) * t123 - g(3) * t120;
t144 = -t111 * t29 - t81;
t36 = pkin(4) * t179 + t50;
t41 = t47 - t206;
t143 = t110 * t41 + t111 * t36;
t124 = qJD(4) ^ 2;
t90 = -pkin(3) - t207;
t142 = t110 * t90 + t124 * t89;
t17 = t177 + t139;
t141 = t13 * t118 - t17 * t121;
t140 = t23 * t118 - t24 * t121;
t138 = -t147 - t8;
t137 = -t14 - t147;
t136 = g(1) * t118 + t121 * t82 - t5;
t135 = t146 + t97;
t134 = -qJDD(4) * t89 + t90 * t180;
t132 = t147 + t172;
t131 = t110 * t47 + t124 * t48 + t193;
t130 = -qJDD(4) * t48 + (t111 * t47 - t52) * qJD(4);
t129 = -t81 - t189 + (-qJD(5) - t22) * t111;
t126 = t201 + t127;
t107 = t123 * pkin(1);
t106 = t120 * pkin(1);
t105 = t121 * qJ(5);
t72 = t118 * t187;
t64 = qJDD(4) * t121 - t124 * t118;
t63 = qJDD(4) * t118 + t124 * t121;
t58 = t182 * t109;
t57 = t121 * t89 + t105;
t56 = t190 * t118;
t46 = t113 * t110 - 0.2e1 * t150;
t45 = t112 * t110 + 0.2e1 * t150;
t43 = -t118 * qJD(5) + t121 * t158;
t42 = t118 * t158 + t101;
t38 = t51 * t178;
t32 = t121 * t48 + t105;
t31 = t191 * t118;
t28 = 0.2e1 * t118 * t86 - 0.2e1 * t182 * t180;
t25 = t29 * t179;
t18 = t22 * t179;
t10 = (-qJD(5) - t52) * t118 + t121 * t159;
t9 = t118 * t159 + t121 * t52 + t101;
t2 = t188 + t99 - qJD(4) * t139 + (-qJD(5) * t111 - t189 - t213) * t118;
t1 = [0, 0, 0, 0, 0, qJDD(1), t145, g(2) * t120 - g(3) * t123, 0, 0, 0, 0, 0, 0, 0, t110, (t110 * t122 + t119 * t154) * pkin(1) + t135, ((-qJDD(1) - t110) * t119 + t122 * t154) * pkin(1) + t164, 0, (t145 + (t119 ^ 2 + t122 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t110, t53 * t110 - t132 - t193, -t54 * t110 + t148 - t192, 0, t21 * t54 + t35 * t52 - t172 * t53 - t34 * t50 - g(2) * (t94 + t107) - g(3) * (t93 + t106), t45, t28, t63, t46, t64, 0, t25 + t130 * t118 + (-t131 + t137) * t121, t131 * t118 + t130 * t121 + t161, t48 * t157 + t181 * t192 + t126, t14 * t47 + t29 * t50 - g(2) * (t107 + t171) - g(3) * (t106 + t153) - t140 * t52 + t127 * t48, t45, t28, t63, t46, t64, 0, t31 * qJDD(4) + t18 + (t41 * t186 + t10) * qJD(4) + (t138 - t143) * t121, -t32 * qJDD(4) + t143 * t118 + (t41 * t185 - t9) * qJD(4) + t163, (t110 * t32 + t111 * t9 + (-t111 * t31 - t13) * qJD(4)) * t121 + (-t10 * t111 - t110 * t31 - t2 + (-t111 * t32 - t17) * qJD(4)) * t118 + t174, t3 * t32 + t17 * t9 + t2 * t31 + t13 * t10 + t8 * t41 + t22 * t36 - g(2) * (t107 + t149) - g(3) * (t106 + t173); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t119 * pkin(1) * t155 + t135, (t122 * t155 - t175) * pkin(1) + t164, 0, 0, 0, 0, 0, 0, 0, t110, t110 * t207 - t132 + t194, -t110 * t208 + t51 * t111 + t148, 0, t34 * t49 - t35 * t51 + (t115 * t21 - t116 * t172 + t146) * pkin(2), t45, t28, t63, t46, t64, 0, t25 + t134 * t118 + (t137 - t142) * t121 + t203, t38 + t134 * t121 + (t142 - t194) * t118 + t161, -t111 * t162 + t89 * t157 + t126, -g(2) * t171 - g(3) * t153 + t127 * t89 + t14 * t90 + t140 * t51 - t29 * t49, t45, t28, t63, t46, t64, 0, t56 * qJDD(4) + t18 + (t62 * t186 + t43) * qJD(4) + (t138 - t74 - t197) * t121 + t203, -t57 * qJDD(4) + t38 + (-t194 + t197) * t118 + (-t42 + (t121 * t62 + t210) * t111) * qJD(4) + t163, (-qJD(4) * t13 + t110 * t57) * t121 + (-t17 * qJD(4) - t110 * t56 - t2) * t118 + (-t118 * t43 + t121 * t42 - t162 + (-t118 * t57 - t121 * t56) * qJD(4)) * t111 + t174, t3 * t57 + t2 * t56 + t13 * t43 + t8 * t62 - t22 * t49 - g(2) * t149 - g(3) * t173 + (-t121 * t51 + t42) * t17 + (t13 * t51 + t22 * t198) * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, 0, 0, 0, t64, -t63, 0, -t140 * qJD(4) + t5 * t118 + t6 * t121 - g(1), 0, 0, 0, 0, 0, 0, t64, -t63, 0, -qJD(4) * t141 + t3 * t118 + t2 * t121 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t58, t85, t72, t86, qJDD(4), (t24 - t195) * qJD(4) + (t144 - t213) * t118 + t151, t23 * qJD(4) + t144 * t121 + t136, 0, 0, -t72, t58, t85, t72, t86, qJDD(4), 0.2e1 * t188 + (t17 - t139) * qJD(4) + (pkin(4) * t187 + t129 - t213) * t118 + t151, -t109 * t210 + (qJ(5) * t186 + t16) * qJD(4) + t129 * t121 + t136, -pkin(4) * t85 + (-t198 + t204) * t185, t204 * t17 + (-t209 + t2 + (-t111 * t22 - t201) * t118) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86 + 0.2e1 * t166, t85 + 0.2e1 * t165, -t181 * t109, t111 * t141 - t138;];
tau_reg = t1;
