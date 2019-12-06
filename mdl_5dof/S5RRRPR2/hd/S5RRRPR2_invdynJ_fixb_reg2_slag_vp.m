% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:41:04
% EndTime: 2019-12-05 18:41:06
% DurationCPUTime: 1.45s
% Computational Cost: add. (3105->246), mult. (5261->303), div. (0->0), fcn. (3029->16), ass. (0->161)
t107 = sin(qJ(3));
t111 = cos(qJ(3));
t108 = sin(qJ(2));
t159 = qJDD(1) * t108;
t112 = cos(qJ(2));
t163 = qJD(2) * t112;
t127 = (qJD(1) * t163 + t159) * pkin(1);
t179 = pkin(1) * qJD(1);
t156 = t108 * t179;
t189 = pkin(1) * t112;
t91 = qJDD(1) * t189;
t99 = qJDD(1) + qJDD(2);
t49 = pkin(2) * t99 - qJD(2) * t156 + t91;
t100 = qJD(1) + qJD(2);
t64 = pkin(2) * t100 + t112 * t179;
t204 = -t107 * t49 - (qJD(3) * t64 + t127) * t111;
t103 = qJ(1) + qJ(2);
t98 = qJ(3) + t103;
t86 = pkin(9) + t98;
t75 = sin(t86);
t76 = cos(t86);
t203 = g(2) * t76 + g(3) * t75;
t114 = qJD(5) ^ 2;
t104 = sin(pkin(9));
t105 = cos(pkin(9));
t170 = t105 * t107;
t178 = pkin(2) * qJD(3);
t167 = t108 * t111;
t135 = -t107 * t112 - t167;
t57 = t135 * t179;
t168 = t107 * t108;
t134 = t111 * t112 - t168;
t58 = t134 * t179;
t182 = -t104 * t58 + t105 * t57 + (t104 * t111 + t170) * t178;
t94 = qJD(3) + t100;
t153 = t182 * t94;
t171 = t104 * t107;
t89 = pkin(2) * t111 + pkin(3);
t55 = -pkin(2) * t171 + t105 * t89;
t50 = -pkin(4) - t55;
t56 = pkin(2) * t170 + t104 * t89;
t51 = pkin(8) + t56;
t93 = qJDD(3) + t99;
t202 = t114 * t51 + t50 * t93 + t153;
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t41 = -t107 * t156 + t111 * t64;
t38 = pkin(3) * t94 + t41;
t42 = t107 * t64 + t111 * t156;
t39 = t105 * t42;
t24 = t104 * t38 + t39;
t20 = pkin(8) * t94 + t24;
t13 = qJD(4) * t110 - t106 * t20;
t175 = t110 * t20;
t14 = qJD(4) * t106 + t175;
t173 = qJD(5) * t13;
t161 = qJD(3) * t111;
t152 = t108 * t161;
t162 = qJD(3) * t107;
t43 = t111 * t49;
t22 = -t64 * t162 + t43 + (-t107 * t159 + (-t107 * t163 - t152) * qJD(1)) * pkin(1);
t15 = pkin(3) * t93 + t22;
t145 = qJD(3) * t156;
t67 = t107 * t145;
t21 = -t204 - t67;
t8 = t104 * t15 + t105 * t21;
t6 = pkin(8) * t93 + t8;
t2 = qJDD(4) * t106 + t110 * t6 + t173;
t172 = qJD(5) * t14;
t95 = t110 * qJDD(4);
t3 = -t106 * t6 - t172 + t95;
t120 = -t3 * t106 + t2 * t110 + (-t106 * t14 - t110 * t13) * qJD(5);
t180 = g(2) * t75 - g(3) * t76;
t87 = sin(t98);
t88 = cos(t98);
t200 = g(2) * t88 + g(3) * t87;
t96 = sin(t103);
t97 = cos(t103);
t199 = g(2) * t97 + g(3) * t96;
t181 = -t104 * t57 - t105 * t58 + (t105 * t111 - t171) * t178;
t198 = t181 * t94;
t101 = t106 ^ 2;
t102 = t110 ^ 2;
t164 = t101 + t102;
t197 = pkin(2) * t93;
t196 = pkin(2) * t96;
t195 = pkin(2) * t97;
t194 = pkin(3) * t87;
t193 = pkin(3) * t88;
t90 = pkin(2) + t189;
t34 = t90 * t161 + (t134 * qJD(2) - t108 * t162) * pkin(1);
t35 = -t90 * t162 + (t135 * qJD(2) - t152) * pkin(1);
t9 = t104 * t34 - t105 * t35;
t192 = t9 * t94;
t160 = qJD(5) * t110;
t177 = t104 * t42;
t23 = t105 * t38 - t177;
t19 = -pkin(4) * t94 - t23;
t151 = t104 * t21 - t105 * t15;
t5 = -pkin(4) * t93 + t151;
t191 = t5 * t106 + t19 * t160;
t109 = sin(qJ(1));
t190 = pkin(1) * t109;
t113 = cos(qJ(1));
t188 = pkin(1) * t113;
t187 = pkin(3) * t104;
t186 = pkin(3) * t105;
t10 = t104 * t35 + t105 * t34;
t185 = t10 * t94;
t25 = t104 * t41 + t39;
t184 = t25 * t94;
t26 = t105 * t41 - t177;
t183 = t26 * t94;
t59 = -pkin(1) * t168 + t111 * t90;
t54 = pkin(3) + t59;
t60 = pkin(1) * t167 + t107 * t90;
t32 = t104 * t54 + t105 * t60;
t169 = t106 * t110;
t166 = qJDD(4) - g(1);
t165 = t101 - t102;
t158 = t19 * qJD(5) * t106 + t203 * t110;
t157 = t91 + t199;
t92 = t94 ^ 2;
t155 = t92 * t169;
t154 = -g(2) * t96 + g(3) * t97;
t150 = t164 * t93;
t149 = -t180 - t8;
t148 = qJD(1) * (-qJD(2) + t100);
t147 = qJD(2) * (-qJD(1) - t100);
t146 = -g(2) * t87 + g(3) * t88 + t67;
t144 = t106 * t94 * t160;
t142 = -t190 - t196;
t141 = -t188 - t195;
t140 = g(2) * t113 + g(3) * t109;
t31 = -t104 * t60 + t105 * t54;
t137 = t106 * t13 - t110 * t14;
t136 = -t151 + t203;
t133 = -pkin(4) * t75 + t76 * pkin(8) - t194;
t132 = -pkin(4) * t76 - pkin(8) * t75 - t193;
t27 = -pkin(4) - t31;
t28 = pkin(8) + t32;
t130 = -t114 * t28 - t27 * t93 - t192;
t81 = pkin(8) + t187;
t82 = -pkin(4) - t186;
t129 = -t114 * t81 - t82 * t93 + t184;
t128 = t133 - t196;
t126 = -qJDD(5) * t28 + (t27 * t94 - t10) * qJD(5);
t125 = -qJDD(5) * t81 + (t82 * t94 + t26) * qJD(5);
t124 = t132 - t195;
t123 = -qJDD(5) * t51 + (t50 * t94 - t181) * qJD(5);
t121 = -qJD(4) * qJD(5) - t19 * t94 - t180 - t6;
t119 = t180 + t120;
t118 = (-pkin(2) * t94 - t64) * qJD(3) - t127;
t117 = t146 + t204;
t116 = t22 + t200;
t66 = qJDD(5) * t110 - t106 * t114;
t65 = qJDD(5) * t106 + t110 * t114;
t45 = t102 * t93 - 0.2e1 * t144;
t44 = t101 * t93 + 0.2e1 * t144;
t36 = -0.2e1 * t165 * t94 * qJD(5) + 0.2e1 * t93 * t169;
t1 = [0, 0, 0, 0, 0, qJDD(1), t140, -g(2) * t109 + g(3) * t113, 0, 0, 0, 0, 0, 0, 0, t99, (t108 * t147 + t112 * t99) * pkin(1) + t157, ((-qJDD(1) - t99) * t108 + t112 * t147) * pkin(1) + t154, 0, (t140 + (t108 ^ 2 + t112 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t93, t35 * t94 + t59 * t93 + t116, -t34 * t94 - t60 * t93 + t117, 0, -g(2) * t141 - g(3) * t142 + t21 * t60 + t22 * t59 + t42 * t34 + t41 * t35, 0, 0, 0, 0, 0, t93, t31 * t93 + t136 - t192, -t32 * t93 + t149 - t185, 0, t8 * t32 + t24 * t10 - t151 * t31 - t23 * t9 - g(2) * (t141 - t193) - g(3) * (t142 - t194), t44, t36, t65, t45, t66, 0, t126 * t106 + (t130 - t5) * t110 + t158, t126 * t110 + (-t130 - t203) * t106 + t191, t28 * t150 + t164 * t185 + t119, t5 * t27 + t19 * t9 - g(2) * (t124 - t188) - g(3) * (t128 - t190) - t137 * t10 + t120 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t108 * pkin(1) * t148 + t157, (t112 * t148 - t159) * pkin(1) + t154, 0, 0, 0, 0, 0, 0, 0, t93, -t57 * t94 + t43 + (-t145 + t197) * t111 + t118 * t107 + t200, t58 * t94 + (-t49 - t197) * t107 + t118 * t111 + t146, 0, -t41 * t57 - t42 * t58 + (t107 * t21 + t111 * t22 + (-t107 * t41 + t111 * t42) * qJD(3) + t199) * pkin(2), 0, 0, 0, 0, 0, t93, t55 * t93 + t136 - t153, -t56 * t93 + t149 - t198, 0, t8 * t56 - t151 * t55 - g(2) * (-t193 - t195) - g(3) * (-t194 - t196) + t181 * t24 - t182 * t23, t44, t36, t65, t45, t66, 0, t123 * t106 + (-t5 - t202) * t110 + t158, t123 * t110 + (-t203 + t202) * t106 + t191, t51 * t150 + t164 * t198 + t119, t5 * t50 - g(2) * t124 - g(3) * t128 + t182 * t19 + ((t2 - t173) * t51 + t181 * t14) * t110 + ((-t3 - t172) * t51 - t181 * t13) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t42 * t94 + t116, t41 * t94 + t117, 0, 0, 0, 0, 0, 0, 0, t93, t93 * t186 + t136 + t184, -t93 * t187 + t149 + t183, 0, t23 * t25 - t24 * t26 + (t104 * t8 - t105 * t151 + t200) * pkin(3), t44, t36, t65, t45, t66, 0, t125 * t106 + (t129 - t5) * t110 + t158, t125 * t110 + (-t129 - t203) * t106 + t191, t119 + t164 * (t81 * t93 - t183), -g(2) * t132 - g(3) * t133 + t120 * t81 + t137 * t26 - t19 * t25 + t5 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, 0, 0, 0, 0, 0, 0, t66, -t65, 0, -qJD(5) * t137 + t106 * t2 + t110 * t3 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, t165 * t92, t106 * t93, t155, t110 * t93, qJDD(5), -g(1) * t110 + t95 + (t14 - t175) * qJD(5) + t121 * t106, t173 + (qJD(5) * t20 - t166) * t106 + t121 * t110, 0, 0;];
tau_reg = t1;
