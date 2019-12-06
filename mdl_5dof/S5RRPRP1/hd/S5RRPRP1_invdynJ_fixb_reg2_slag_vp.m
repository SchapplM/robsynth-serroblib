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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:22:30
% EndTime: 2019-12-05 18:22:33
% DurationCPUTime: 1.60s
% Computational Cost: add. (1918->281), mult. (3139->325), div. (0->0), fcn. (1791->12), ass. (0->171)
t106 = qJDD(1) + qJDD(2);
t111 = sin(pkin(8));
t112 = cos(pkin(8));
t118 = cos(qJ(2));
t192 = pkin(1) * qJD(2);
t162 = qJD(1) * t192;
t115 = sin(qJ(2));
t168 = qJDD(1) * t115;
t210 = pkin(1) * t168 + t118 * t162;
t206 = pkin(1) * t118;
t95 = qJDD(1) * t206;
t44 = pkin(2) * t106 - t115 * t162 + t95;
t21 = t111 * t44 + t210 * t112;
t15 = pkin(7) * t106 + t21;
t211 = qJD(3) * qJD(4) + t15;
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t100 = t117 * qJD(3);
t107 = qJD(1) + qJD(2);
t193 = pkin(1) * qJD(1);
t164 = t115 * t193;
t163 = t118 * t193;
t59 = pkin(2) * t107 + t163;
t35 = t111 * t59 + t112 * t164;
t30 = pkin(7) * t107 + t35;
t23 = -t114 * t30 + t100;
t170 = t114 * qJD(3);
t186 = t117 * t30;
t24 = t170 + t186;
t172 = qJD(4) * t114;
t5 = t114 * qJDD(3) + t211 * t117 - t30 * t172;
t97 = t117 * qJDD(3);
t6 = -t24 * qJD(4) - t114 * t15 + t97;
t123 = -t6 * t114 + t5 * t117 + (-t114 * t24 - t117 * t23) * qJD(4);
t110 = qJ(1) + qJ(2);
t98 = pkin(8) + t110;
t85 = sin(t98);
t77 = g(3) * t85;
t86 = cos(t98);
t143 = -g(2) * t86 - t77;
t78 = g(3) * t86;
t79 = g(2) * t85;
t194 = t79 - t78;
t101 = sin(t110);
t102 = cos(t110);
t209 = g(2) * t102 + g(3) * t101;
t155 = qJ(5) * t107 + t30;
t136 = t155 * t117;
t171 = qJD(4) * t117;
t75 = t111 * t164;
t34 = t112 * t59 - t75;
t199 = pkin(4) * t117;
t93 = pkin(3) + t199;
t22 = -t93 * t107 + qJD(5) - t34;
t160 = t107 * t172;
t166 = t210 * t111 - t112 * t44;
t8 = pkin(4) * t160 - t93 * t106 + qJDD(5) + t166;
t208 = t8 * t114 + t22 * t171;
t116 = sin(qJ(1));
t207 = pkin(1) * t116;
t119 = cos(qJ(1));
t205 = pkin(1) * t119;
t204 = pkin(2) * t101;
t203 = pkin(2) * t102;
t202 = pkin(2) * t111;
t201 = pkin(2) * t112;
t108 = t114 ^ 2;
t200 = pkin(4) * t108;
t198 = g(1) * t117;
t14 = -t106 * pkin(3) + t166;
t29 = -pkin(3) * t107 - t34;
t197 = t14 * t114 + t29 * t171;
t16 = -t155 * t114 + t100;
t191 = qJD(4) * pkin(4);
t13 = t16 + t191;
t196 = t13 - t16;
t185 = t117 * t86;
t195 = g(2) * t185 + t117 * t77;
t176 = t112 * t115;
t94 = pkin(2) + t206;
t54 = pkin(1) * t176 + t111 * t94;
t60 = -t93 - t201;
t190 = t106 * t60;
t131 = pkin(1) * (t111 * t118 + t176);
t49 = qJD(1) * t131;
t189 = t107 * t49;
t50 = qJD(2) * t131;
t188 = t107 * t50;
t177 = t111 * t115;
t52 = (t112 * t118 - t177) * t192;
t187 = t107 * t52;
t48 = pkin(7) + t54;
t184 = -qJ(5) - t48;
t87 = pkin(7) + t202;
t183 = -qJ(5) - t87;
t182 = qJ(5) * t106;
t181 = qJDD(4) * pkin(4);
t105 = t107 ^ 2;
t180 = t105 * t117;
t179 = t107 * t114;
t178 = t107 * t117;
t83 = t114 * t106;
t84 = t117 * t106;
t109 = t117 ^ 2;
t175 = t108 - t109;
t174 = t108 + t109;
t173 = qJD(4) * t107;
t99 = t117 * qJD(5);
t3 = t107 * t99 + (-t160 + t84) * qJ(5) + t5;
t167 = t3 * t117 + t194;
t165 = t95 + t209;
t159 = t107 * t171;
t158 = -g(2) * t101 + g(3) * t102;
t53 = -pkin(1) * t177 + t112 * t94;
t51 = t112 * t163 - t75;
t157 = t174 * t51;
t156 = t51 * t172 + t49 * t178 + t195;
t154 = qJD(4) * t184;
t153 = qJD(4) * t183;
t152 = t106 * t174;
t150 = qJD(1) * (-qJD(2) + t107);
t149 = qJD(2) * (-qJD(1) - t107);
t47 = -pkin(3) - t53;
t147 = t114 * t78 - t198 + t97;
t146 = t114 * t159;
t145 = t166 + t143;
t144 = -t21 - t194;
t142 = g(2) * t119 + g(3) * t116;
t141 = -t107 * t29 - t79;
t36 = pkin(4) * t172 + t50;
t41 = t47 - t199;
t140 = -t106 * t41 - t107 * t36;
t120 = qJD(4) ^ 2;
t88 = -pkin(3) - t201;
t139 = -t106 * t88 - t120 * t87;
t17 = t170 + t136;
t138 = t114 * t13 - t117 * t17;
t137 = t114 * t23 - t117 * t24;
t135 = g(1) * t114 + g(3) * t185 - t5;
t134 = -pkin(3) * t85 + t86 * pkin(7) - t204;
t113 = -qJ(5) - pkin(7);
t133 = t85 * t113 - t86 * t93 - t203;
t132 = -qJDD(4) * t87 + t88 * t173;
t130 = -pkin(3) * t86 - pkin(7) * t85 - t203;
t129 = t143 - t189;
t128 = -t113 * t86 - t85 * t93 - t204;
t127 = -t106 * t47 - t120 * t48 - t188;
t126 = -qJDD(4) * t48 + (t107 * t47 - t52) * qJD(4);
t125 = -t79 - t182 + (-qJD(5) - t22) * t107;
t122 = t194 + t123;
t103 = t117 * qJ(5);
t70 = t114 * t180;
t62 = qJDD(4) * t117 - t114 * t120;
t61 = qJDD(4) * t114 + t117 * t120;
t58 = t175 * t105;
t57 = t117 * t87 + t103;
t56 = t183 * t114;
t46 = t106 * t109 - 0.2e1 * t146;
t45 = t106 * t108 + 0.2e1 * t146;
t43 = -t114 * qJD(5) + t117 * t153;
t42 = t114 * t153 + t99;
t38 = t51 * t171;
t32 = t117 * t48 + t103;
t31 = t184 * t114;
t28 = 0.2e1 * t114 * t84 - 0.2e1 * t175 * t173;
t25 = t29 * t172;
t18 = t22 * t172;
t10 = (-qJD(5) - t52) * t114 + t117 * t154;
t9 = t114 * t154 + t117 * t52 + t99;
t2 = t181 + t97 - qJD(4) * t136 + (-qJD(5) * t107 - t182 - t211) * t114;
t1 = [0, 0, 0, 0, 0, qJDD(1), t142, -g(2) * t116 + g(3) * t119, 0, 0, 0, 0, 0, 0, 0, t106, (t106 * t118 + t115 * t149) * pkin(1) + t165, ((-qJDD(1) - t106) * t115 + t118 * t149) * pkin(1) + t158, 0, (t142 + (t115 ^ 2 + t118 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t106, t106 * t53 - t145 - t188, -t106 * t54 + t144 - t187, 0, t21 * t54 + t35 * t52 - t166 * t53 - t34 * t50 - g(2) * (-t203 - t205) - g(3) * (-t204 - t207), t45, t28, t61, t46, t62, 0, t25 + t126 * t114 + (t127 - t14) * t117 + t195, t126 * t117 + (-t127 + t143) * t114 + t197, t48 * t152 + t174 * t187 + t122, t14 * t47 + t29 * t50 - g(2) * (t130 - t205) - g(3) * (t134 - t207) - t137 * t52 + t123 * t48, t45, t28, t61, t46, t62, 0, t31 * qJDD(4) + t18 + (t41 * t179 + t10) * qJD(4) + (t140 - t8) * t117 + t195, -t32 * qJDD(4) + (t41 * t178 - t9) * qJD(4) + (-t140 + t143) * t114 + t208, (t106 * t32 + t107 * t9 + (-t107 * t31 - t13) * qJD(4)) * t117 + (-t10 * t107 - t106 * t31 - t2 + (-t107 * t32 - t17) * qJD(4)) * t114 + t167, t3 * t32 + t17 * t9 + t2 * t31 + t13 * t10 + t8 * t41 + t22 * t36 - g(2) * (t133 - t205) - g(3) * (t128 - t207); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t115 * pkin(1) * t150 + t165, (t118 * t150 - t168) * pkin(1) + t158, 0, 0, 0, 0, 0, 0, 0, t106, t106 * t201 - t145 + t189, -t106 * t202 + t107 * t51 + t144, 0, t34 * t49 - t35 * t51 + (t111 * t21 - t112 * t166 + t209) * pkin(2), t45, t28, t61, t46, t62, 0, t25 + t132 * t114 + (t139 - t14) * t117 + t156, t38 + t132 * t117 + (t129 - t139) * t114 + t197, -t107 * t157 + t87 * t152 + t122, -g(2) * t130 - g(3) * t134 + t123 * t87 + t137 * t51 + t14 * t88 - t29 * t49, t45, t28, t61, t46, t62, 0, t56 * qJDD(4) + t18 + (-t8 - t190) * t117 + (t43 + (t60 - t199) * t179) * qJD(4) + t156, -t57 * qJDD(4) + t38 + (-t42 + (t117 * t60 + t200) * t107) * qJD(4) + (t129 + t190) * t114 + t208, (-qJD(4) * t13 + t106 * t57) * t117 + (-t17 * qJD(4) - t106 * t56 - t2) * t114 + (-t114 * t43 + t117 * t42 - t157 + (-t114 * t57 - t117 * t56) * qJD(4)) * t107 + t167, t3 * t57 + t2 * t56 + t13 * t43 + t8 * t60 - t22 * t49 - g(2) * t133 - g(3) * t128 + (-t117 * t51 + t42) * t17 + (t13 * t51 + t22 * t191) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, 0, 0, 0, t62, -t61, 0, -t137 * qJD(4) + t5 * t114 + t6 * t117 - g(1), 0, 0, 0, 0, 0, 0, t62, -t61, 0, -t138 * qJD(4) + t3 * t114 + t2 * t117 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t58, t83, t70, t84, qJDD(4), (t24 - t186) * qJD(4) + (t141 - t211) * t114 + t147, t23 * qJD(4) + t141 * t117 + t135, 0, 0, -t70, t58, t83, t70, t84, qJDD(4), 0.2e1 * t181 + (t17 - t136) * qJD(4) + (pkin(4) * t180 + t125 - t211) * t114 + t147, -t105 * t200 + (qJ(5) * t179 + t16) * qJD(4) + t125 * t117 + t135, -pkin(4) * t83 + (-t191 + t196) * t178, t196 * t17 + (-t198 + t2 + (-t107 * t22 - t194) * t114) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84 + 0.2e1 * t160, t83 + 0.2e1 * t159, -t174 * t105, t138 * t107 + t143 + t8;];
tau_reg = t1;
