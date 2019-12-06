% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:56
% EndTime: 2019-12-05 15:58:02
% DurationCPUTime: 1.71s
% Computational Cost: add. (1525->263), mult. (3715->375), div. (0->0), fcn. (3162->14), ass. (0->151)
t103 = sin(qJ(4));
t97 = sin(pkin(10));
t179 = t103 * t97;
t150 = qJD(2) * t179;
t106 = cos(qJ(4));
t100 = cos(pkin(10));
t166 = qJD(2) * t100;
t85 = t106 * t166;
t62 = -t85 + t150;
t55 = qJD(5) + t62;
t197 = t55 - qJD(5);
t101 = cos(pkin(5));
t107 = cos(qJ(2));
t173 = cos(pkin(9));
t140 = t173 * t107;
t104 = sin(qJ(2));
t98 = sin(pkin(9));
t178 = t104 * t98;
t58 = -t101 * t140 + t178;
t141 = t173 * t104;
t174 = t98 * t107;
t60 = t101 * t174 + t141;
t138 = g(1) * t60 + g(2) * t58;
t99 = sin(pkin(5));
t176 = t107 * t99;
t114 = g(3) * t176 - t138;
t96 = pkin(10) + qJ(4);
t90 = cos(t96);
t111 = t114 * t90;
t157 = t100 * qJDD(2);
t161 = t97 * qJDD(2);
t135 = t103 * t161 - t106 * t157;
t69 = t100 * t103 + t106 * t97;
t66 = t69 * qJD(4);
t32 = qJD(2) * t66 + t135;
t27 = qJDD(5) + t32;
t127 = t106 * t100 - t179;
t87 = -pkin(3) * t100 - pkin(2);
t30 = -pkin(4) * t127 - pkin(8) * t69 + t87;
t196 = t30 * t27 - t111;
t108 = qJD(2) ^ 2;
t119 = (qJDD(2) * t107 - t104 * t108) * t99;
t64 = t69 * qJD(2);
t172 = qJD(1) * t99;
t151 = t107 * t172;
t133 = qJD(3) - t151;
t152 = t104 * t172;
t71 = qJD(2) * qJ(3) + t152;
t167 = qJD(1) * t101;
t81 = t100 * t167;
t37 = t81 + (-pkin(7) * qJD(2) - t71) * t97;
t50 = t100 * t71 + t97 * t167;
t38 = pkin(7) * t166 + t50;
t128 = t103 * t38 - t106 * t37;
t158 = qJDD(2) * qJ(3);
t163 = qJDD(1) * t99;
t48 = t104 * t163 + t158 + (qJD(3) + t151) * qJD(2);
t159 = qJDD(1) * t101;
t79 = t100 * t159;
t18 = t79 + (-pkin(7) * qJDD(2) - t48) * t97;
t25 = t100 * t48 + t97 * t159;
t19 = pkin(7) * t157 + t25;
t129 = t103 * t18 + t106 * t19;
t1 = qJDD(4) * pkin(8) - qJD(4) * t128 + t129;
t59 = t101 * t141 + t174;
t61 = -t101 * t178 + t140;
t137 = g(1) * t61 + g(2) * t59;
t54 = qJD(2) * t87 + t133;
t15 = pkin(4) * t62 - pkin(8) * t64 + t54;
t177 = t104 * t99;
t156 = g(3) * t177;
t116 = t127 * t176;
t184 = pkin(7) + qJ(3);
t72 = t184 * t97;
t73 = t184 * t100;
t41 = t103 * t73 + t106 * t72;
t183 = qJD(1) * t116 - qJD(3) * t127 + qJD(4) * t41;
t130 = t103 * t19 - t106 * t18;
t14 = t103 * t37 + t106 * t38;
t2 = -qJDD(4) * pkin(4) + qJD(4) * t14 + t130;
t42 = -t103 * t72 + t106 * t73;
t65 = t127 * qJD(4);
t9 = -qJD(4) * pkin(4) + t128;
t195 = (qJD(5) * t15 + t1) * t127 + t2 * t69 + t9 * t65 + (-qJD(5) * t30 + t183) * t55 - t42 * t27 - t156 - t137;
t170 = qJDD(2) * pkin(2);
t118 = t138 + t170;
t165 = qJD(2) * t104;
t146 = qJD(1) * t165;
t160 = t99 * t146 + qJDD(3);
t122 = -t107 * t163 + t160;
t51 = t122 - t170;
t194 = (-g(3) * t107 + t146) * t99 + t118 - t51;
t147 = t99 * t173;
t185 = t98 * t99;
t89 = sin(t96);
t120 = g(1) * (t90 * t185 - t61 * t89) + g(2) * (-t147 * t90 - t59 * t89) + g(3) * (t101 * t90 - t89 * t177);
t193 = (pkin(4) * t64 + t55 * pkin(8)) * t55 + t120 + t2;
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t155 = qJD(4) * t85 + t103 * t157 + t106 * t161;
t31 = -qJD(4) * t150 + t155;
t47 = qJD(4) * t102 + t105 * t64;
t8 = qJD(5) * t47 - t105 * qJDD(4) + t102 * t31;
t132 = t100 * t50 - (-t71 * t97 + t81) * t97;
t192 = t107 * t132 + (-qJD(2) * pkin(2) + t133) * t104;
t162 = t105 * qJD(4);
t164 = qJD(5) * t102;
t7 = qJD(5) * t162 + t102 * qJDD(4) + t105 * t31 - t164 * t64;
t191 = t102 * t7;
t45 = t102 * t64 - t162;
t189 = t45 * t55;
t188 = t45 * t64;
t187 = t47 * t55;
t186 = t47 * t64;
t115 = t69 * t176;
t182 = -qJD(1) * t115 + qJD(3) * t69 + qJD(4) * t42;
t181 = t100 ^ 2 + t97 ^ 2;
t180 = t102 * t27;
t139 = t105 * t55;
t171 = qJD(5) * t69;
t169 = t107 * t108;
t168 = qJDD(1) - g(3);
t154 = t102 * t176;
t153 = t105 * t176;
t149 = t99 * t165;
t148 = t69 * t164;
t142 = t168 * t99;
t136 = pkin(4) * t66 - pkin(8) * t65 - t152;
t134 = t27 * t69 + t55 * t65;
t10 = qJD(4) * pkin(8) + t14;
t4 = t10 * t105 + t102 * t15;
t131 = t10 * t102 - t105 * t15;
t56 = t100 * t101 - t97 * t177;
t57 = t100 * t177 + t101 * t97;
t20 = t103 * t57 - t106 * t56;
t21 = t103 * t56 + t106 * t57;
t126 = t105 * t27 + (-t102 * t62 - t164) * t55;
t124 = -t102 * t21 - t153;
t123 = -t105 * t21 + t154;
t112 = -pkin(8) * t27 + (-t128 + t9) * t55;
t24 = -t48 * t97 + t79;
t110 = t100 * t25 - t24 * t97 - t137;
t39 = qJDD(2) * t87 + t122;
t53 = t101 * t89 + t90 * t177;
t36 = t89 * t185 + t61 * t90;
t34 = -t147 * t89 + t59 * t90;
t12 = qJD(2) * t115 + qJD(4) * t21;
t11 = qJD(2) * t116 - qJD(4) * t20;
t6 = pkin(4) * t32 - pkin(8) * t31 + t39;
t5 = t105 * t6;
t3 = [t168, 0, t119, (-qJDD(2) * t104 - t169) * t99, t100 * t119, -t97 * t119, t181 * t99 * t169 + (t100 * t57 - t56 * t97) * qJDD(2), t24 * t56 + t25 * t57 - g(3) + (t192 * qJD(2) - t107 * t51) * t99, 0, 0, 0, 0, 0, -qJD(4) * t12 - qJDD(4) * t20 + (-t107 * t32 + t165 * t62) * t99, -qJD(4) * t11 - qJDD(4) * t21 + (-t107 * t31 + t165 * t64) * t99, 0, 0, 0, 0, 0, (qJD(5) * t123 - t102 * t11 + t105 * t149) * t55 + t124 * t27 + t12 * t45 + t20 * t8, -(qJD(5) * t124 + t102 * t149 + t105 * t11) * t55 + t123 * t27 + t12 * t47 + t20 * t7; 0, qJDD(2), t168 * t176 + t138, -t104 * t142 + t137, t194 * t100, -t194 * t97, -t156 + t110 + (t133 * qJD(2) + t158) * t181, t132 * qJD(3) + (t138 - t51) * pkin(2) + t110 * qJ(3) + (-g(3) * (pkin(2) * t107 + qJ(3) * t104) - t192 * qJD(1)) * t99, t31 * t69 + t64 * t65, t127 * t31 - t32 * t69 - t62 * t65 - t64 * t66, qJD(4) * t65 + qJDD(4) * t69, -qJD(4) * t66 + qJDD(4) * t127, 0, -t182 * qJD(4) - qJDD(4) * t41 - t127 * t39 - t62 * t152 + t32 * t87 + t54 * t66 - t111, t183 * qJD(4) - qJDD(4) * t42 + t114 * t89 - t64 * t152 + t31 * t87 + t39 * t69 + t54 * t65, -t47 * t148 + (t47 * t65 + t69 * t7) * t105, (-t102 * t47 - t105 * t45) * t65 + (-t191 - t105 * t8 + (t102 * t45 - t105 * t47) * qJD(5)) * t69, t105 * t134 - t127 * t7 - t148 * t55 + t47 * t66, -t102 * t134 + t127 * t8 - t171 * t139 - t45 * t66, -t127 * t27 + t55 * t66, -t131 * t66 + t41 * t8 - t5 * t127 + t182 * t45 + (t136 * t55 + (t10 * t127 - t42 * t55 + t69 * t9) * qJD(5) + t196) * t105 + t195 * t102, -t4 * t66 + t41 * t7 + t182 * t47 + ((-qJD(5) * t10 + t6) * t127 - t9 * t171 + (qJD(5) * t42 - t136) * t55 - t196) * t102 + t195 * t105; 0, 0, 0, 0, -t157, t161, -t181 * t108, -qJD(2) * t132 - t107 * t142 - t118 + t160, 0, 0, 0, 0, 0, 0.2e1 * qJD(4) * t64 + t135, (-t62 - t150) * qJD(4) + t155, 0, 0, 0, 0, 0, t126 - t188, -t55 ^ 2 * t105 - t180 - t186; 0, 0, 0, 0, 0, 0, 0, 0, t64 * t62, -t62 ^ 2 + t64 ^ 2, (t62 - t150) * qJD(4) + t155, -t135, qJDD(4), -t54 * t64 - t120 - t130, g(1) * t36 + g(2) * t34 + g(3) * t53 + t54 * t62 - t129, t139 * t47 + t191, (t7 - t189) * t105 + (-t8 - t187) * t102, t139 * t55 + t180 - t186, t126 + t188, -t55 * t64, -pkin(4) * t8 + t112 * t102 - t193 * t105 + t131 * t64 - t14 * t45, -pkin(4) * t7 + t193 * t102 + t112 * t105 - t14 * t47 + t4 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t45, -t45 ^ 2 + t47 ^ 2, t7 + t189, t187 - t8, t27, -t102 * t1 + t5 - t9 * t47 - g(1) * (-t102 * t36 + t105 * t60) - g(2) * (-t102 * t34 + t105 * t58) - g(3) * (-t102 * t53 - t153) + t197 * t4, -t105 * t1 - t102 * t6 + t9 * t45 - g(1) * (-t102 * t60 - t105 * t36) - g(2) * (-t102 * t58 - t105 * t34) - g(3) * (-t105 * t53 + t154) - t197 * t131;];
tau_reg = t3;
