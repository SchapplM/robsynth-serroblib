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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:07:41
% EndTime: 2020-01-03 12:07:43
% DurationCPUTime: 1.37s
% Computational Cost: add. (3105->245), mult. (5261->301), div. (0->0), fcn. (3029->16), ass. (0->163)
t109 = sin(qJ(3));
t113 = cos(qJ(3));
t110 = sin(qJ(2));
t162 = qJDD(1) * t110;
t114 = cos(qJ(2));
t166 = qJD(2) * t114;
t129 = (qJD(1) * t166 + t162) * pkin(1);
t101 = qJDD(1) + qJDD(2);
t183 = pkin(1) * qJD(1);
t160 = t110 * t183;
t193 = pkin(1) * t114;
t91 = qJDD(1) * t193;
t49 = pkin(2) * t101 - qJD(2) * t160 + t91;
t102 = qJD(1) + qJD(2);
t64 = pkin(2) * t102 + t114 * t183;
t201 = -t109 * t49 - (qJD(3) * t64 + t129) * t113;
t105 = qJ(1) + qJ(2);
t98 = qJ(3) + t105;
t86 = pkin(9) + t98;
t75 = sin(t86);
t76 = cos(t86);
t145 = -g(2) * t76 - g(3) * t75;
t106 = sin(pkin(9));
t107 = cos(pkin(9));
t164 = qJD(3) * t113;
t157 = t110 * t164;
t165 = qJD(3) * t109;
t43 = t113 * t49;
t22 = -t64 * t165 + t43 + (-t109 * t162 + (-t109 * t166 - t157) * qJD(1)) * pkin(1);
t93 = qJDD(3) + t101;
t15 = pkin(3) * t93 + t22;
t147 = qJD(3) * t160;
t67 = t109 * t147;
t21 = -t201 - t67;
t155 = t106 * t21 - t107 * t15;
t5 = -pkin(4) * t93 + t155;
t137 = t145 - t5;
t108 = sin(qJ(5));
t112 = cos(qJ(5));
t41 = -t109 * t160 + t113 * t64;
t94 = qJD(3) + t102;
t38 = pkin(3) * t94 + t41;
t42 = t109 * t64 + t113 * t160;
t39 = t107 * t42;
t24 = t106 * t38 + t39;
t20 = pkin(8) * t94 + t24;
t13 = qJD(4) * t112 - t108 * t20;
t179 = t112 * t20;
t14 = qJD(4) * t108 + t179;
t176 = qJD(5) * t13;
t8 = t106 * t15 + t107 * t21;
t6 = pkin(8) * t93 + t8;
t2 = t108 * qJDD(4) + t112 * t6 + t176;
t175 = qJD(5) * t14;
t95 = t112 * qJDD(4);
t3 = -t108 * t6 - t175 + t95;
t122 = -t3 * t108 + t2 * t112 + (-t108 * t14 - t112 * t13) * qJD(5);
t185 = g(2) * t75 - g(3) * t76;
t174 = t106 * t109;
t182 = pkin(2) * qJD(3);
t170 = t110 * t113;
t136 = -t109 * t114 - t170;
t57 = t136 * t183;
t171 = t109 * t110;
t135 = t113 * t114 - t171;
t58 = t135 * t183;
t186 = -t106 * t57 - t107 * t58 + (t107 * t113 - t174) * t182;
t199 = t186 * t94;
t173 = t107 * t109;
t187 = -t106 * t58 + t107 * t57 + (t106 * t113 + t173) * t182;
t198 = t187 * t94;
t103 = t108 ^ 2;
t104 = t112 ^ 2;
t167 = t103 + t104;
t197 = pkin(2) * t93;
t90 = pkin(2) + t193;
t34 = t90 * t164 + (t135 * qJD(2) - t110 * t165) * pkin(1);
t35 = -t90 * t165 + (t136 * qJD(2) - t157) * pkin(1);
t9 = t106 * t34 - t107 * t35;
t194 = t9 * t94;
t192 = pkin(3) * t106;
t191 = pkin(3) * t107;
t10 = t106 * t35 + t107 * t34;
t190 = t10 * t94;
t25 = t106 * t41 + t39;
t189 = t25 * t94;
t181 = t106 * t42;
t26 = t107 * t41 - t181;
t188 = t26 * t94;
t59 = -pkin(1) * t171 + t113 * t90;
t54 = pkin(3) + t59;
t60 = pkin(1) * t170 + t109 * t90;
t32 = t106 * t54 + t107 * t60;
t89 = pkin(2) * t113 + pkin(3);
t56 = pkin(2) * t173 + t106 * t89;
t96 = sin(t105);
t84 = pkin(2) * t96;
t111 = sin(qJ(1));
t99 = t111 * pkin(1);
t184 = t84 + t99;
t115 = cos(qJ(1));
t100 = t115 * pkin(1);
t97 = cos(t105);
t85 = pkin(2) * t97;
t178 = t85 + t100;
t172 = t108 * t112;
t169 = qJDD(4) - g(1);
t168 = t103 - t104;
t163 = qJD(5) * t112;
t88 = cos(t98);
t79 = pkin(3) * t88;
t161 = t76 * pkin(4) + t75 * pkin(8) + t79;
t92 = t94 ^ 2;
t159 = t92 * t172;
t158 = g(2) * t96 - g(3) * t97;
t23 = t107 * t38 - t181;
t19 = -pkin(4) * t94 - t23;
t156 = -t137 * t108 + t19 * t163;
t154 = t167 * t93;
t153 = t85 + t161;
t152 = t185 - t8;
t151 = qJD(1) * (-qJD(2) + t102);
t150 = qJD(2) * (-qJD(1) - t102);
t87 = sin(t98);
t78 = pkin(3) * t87;
t149 = t75 * pkin(4) - pkin(8) * t76 + t78;
t148 = g(2) * t87 - g(3) * t88 + t67;
t146 = t108 * t94 * t163;
t144 = -g(2) * t88 - g(3) * t87;
t143 = -g(2) * t97 - g(3) * t96;
t142 = -g(2) * t115 - g(3) * t111;
t116 = qJD(5) ^ 2;
t55 = -pkin(2) * t174 + t107 * t89;
t50 = -pkin(4) - t55;
t51 = pkin(8) + t56;
t140 = t116 * t51 + t50 * t93;
t139 = t149 + t84;
t31 = -t106 * t60 + t107 * t54;
t138 = t108 * t13 - t112 * t14;
t134 = t143 + t91;
t27 = -pkin(4) - t31;
t28 = pkin(8) + t32;
t132 = t116 * t28 + t27 * t93 + t194;
t81 = pkin(8) + t192;
t82 = -pkin(4) - t191;
t131 = t116 * t81 + t82 * t93 - t189;
t130 = t145 - t155;
t128 = t145 - t198;
t127 = -qJDD(5) * t28 + (t27 * t94 - t10) * qJD(5);
t126 = -qJDD(5) * t81 + (t82 * t94 + t26) * qJD(5);
t125 = -qJDD(5) * t51 + (t50 * t94 - t186) * qJD(5);
t123 = -qJD(4) * qJD(5) - t19 * t94 + t185 - t6;
t121 = -t185 + t122;
t120 = (-pkin(2) * t94 - t64) * qJD(3) - t129;
t119 = t148 + t201;
t118 = t144 + t22;
t66 = qJDD(5) * t112 - t108 * t116;
t65 = qJDD(5) * t108 + t112 * t116;
t45 = t104 * t93 - 0.2e1 * t146;
t44 = t103 * t93 + 0.2e1 * t146;
t36 = -0.2e1 * t168 * t94 * qJD(5) + 0.2e1 * t93 * t172;
t16 = t19 * qJD(5) * t108;
t1 = [0, 0, 0, 0, 0, qJDD(1), t142, g(2) * t111 - g(3) * t115, 0, 0, 0, 0, 0, 0, 0, t101, (t101 * t114 + t110 * t150) * pkin(1) + t134, ((-qJDD(1) - t101) * t110 + t114 * t150) * pkin(1) + t158, 0, (t142 + (t110 ^ 2 + t114 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t93, t35 * t94 + t59 * t93 + t118, -t34 * t94 - t60 * t93 + t119, 0, -g(2) * t178 - g(3) * t184 + t21 * t60 + t22 * t59 + t42 * t34 + t41 * t35, 0, 0, 0, 0, 0, t93, t31 * t93 + t130 - t194, -t32 * t93 + t152 - t190, 0, t8 * t32 + t24 * t10 - t155 * t31 - t23 * t9 - g(2) * (t79 + t178) - g(3) * (t78 + t184), t44, t36, t65, t45, t66, 0, t16 + t127 * t108 + (-t132 + t137) * t112, t132 * t108 + t127 * t112 + t156, t28 * t154 + t167 * t190 + t121, t5 * t27 + t19 * t9 - g(2) * (t100 + t153) - g(3) * (t139 + t99) - t138 * t10 + t122 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t110 * pkin(1) * t151 + t134, (t114 * t151 - t162) * pkin(1) + t158, 0, 0, 0, 0, 0, 0, 0, t93, -t57 * t94 + t43 + (-t147 + t197) * t113 + t120 * t109 + t144, t58 * t94 + (-t49 - t197) * t109 + t120 * t113 + t148, 0, -t41 * t57 - t42 * t58 + (t109 * t21 + t113 * t22 + (-t109 * t41 + t113 * t42) * qJD(3) + t143) * pkin(2), 0, 0, 0, 0, 0, t93, t55 * t93 + t128 - t155, -t56 * t93 + t152 - t199, 0, t8 * t56 - t155 * t55 - g(2) * (t79 + t85) - g(3) * (t78 + t84) + t186 * t24 - t187 * t23, t44, t36, t65, t45, t66, 0, t16 + t125 * t108 + (t128 - t140 - t5) * t112, t125 * t112 + (t140 + t198) * t108 + t156, t51 * t154 + t167 * t199 + t121, t5 * t50 - g(2) * t153 - g(3) * t139 + t187 * t19 + ((t2 - t176) * t51 + t186 * t14) * t112 + ((-t3 - t175) * t51 - t186 * t13) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t42 * t94 + t118, t41 * t94 + t119, 0, 0, 0, 0, 0, 0, 0, t93, t93 * t191 + t130 + t189, -t93 * t192 + t152 + t188, 0, t23 * t25 - t24 * t26 + (t106 * t8 - t107 * t155 + t144) * pkin(3), t44, t36, t65, t45, t66, 0, t16 + t126 * t108 + (-t131 + t137) * t112, t108 * t131 + t112 * t126 + t156, t121 + t167 * (t81 * t93 - t188), -g(2) * t161 - g(3) * t149 + t122 * t81 + t138 * t26 - t19 * t25 + t5 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, 0, 0, 0, 0, 0, 0, t66, -t65, 0, -t138 * qJD(5) + t2 * t108 + t3 * t112 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, t168 * t92, t108 * t93, t159, t112 * t93, qJDD(5), -g(1) * t112 + t95 + (t14 - t179) * qJD(5) + t123 * t108, t176 + (qJD(5) * t20 - t169) * t108 + t123 * t112, 0, 0;];
tau_reg = t1;
