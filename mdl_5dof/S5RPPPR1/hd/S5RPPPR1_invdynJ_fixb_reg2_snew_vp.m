% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:34
% EndTime: 2020-01-03 11:20:44
% DurationCPUTime: 2.04s
% Computational Cost: add. (4582->239), mult. (11021->365), div. (0->0), fcn. (7063->10), ass. (0->158)
t123 = sin(pkin(8));
t117 = t123 ^ 2;
t126 = cos(pkin(8));
t119 = t126 ^ 2;
t132 = qJD(1) ^ 2;
t99 = (t117 + t119) * t132;
t144 = -t126 * pkin(3) - t123 * qJ(4);
t167 = qJD(1) * t123;
t124 = sin(pkin(7));
t127 = cos(pkin(7));
t129 = sin(qJ(1));
t131 = cos(qJ(1));
t147 = -t131 * g(2) - t129 * g(3);
t95 = qJDD(1) * pkin(1) + t147;
t146 = t129 * g(2) - t131 * g(3);
t96 = -t132 * pkin(1) - t146;
t151 = -t124 * t96 + t127 * t95;
t63 = -qJDD(1) * pkin(2) - t132 * qJ(3) + qJDD(3) - t151;
t198 = -0.2e1 * qJD(4) * t167 + t144 * qJDD(1) + t63;
t166 = t126 * qJD(1);
t108 = -qJD(5) + t166;
t122 = sin(pkin(9));
t125 = cos(pkin(9));
t128 = sin(qJ(5));
t130 = cos(qJ(5));
t142 = t122 * t130 + t125 * t128;
t138 = t123 * t142;
t76 = qJD(1) * t138;
t179 = t76 * t108;
t163 = t123 * qJDD(1);
t153 = t122 * t163;
t171 = t123 * t125;
t156 = t130 * t171;
t60 = -t76 * qJD(5) + qJDD(1) * t156 - t128 * t153;
t197 = t60 + t179;
t189 = 2 * qJD(3);
t186 = t124 * t95 + t127 * t96;
t64 = -t132 * pkin(2) + qJDD(1) * qJ(3) + t186;
t196 = qJD(1) * t189 + t64;
t152 = pkin(1) * t124 + qJ(3);
t154 = pkin(1) * t127 + pkin(2);
t195 = -t154 * qJDD(1) + t152 * t99 + t63;
t78 = -t128 * t122 * t167 + qJD(1) * t156;
t194 = (qJD(5) + t108) * t78;
t162 = t126 * qJDD(1);
t107 = -qJDD(5) + t162;
t187 = t78 * t76;
t139 = -t107 - t187;
t193 = t128 * t139;
t192 = t130 * t139;
t74 = t76 ^ 2;
t75 = t78 ^ 2;
t106 = t108 ^ 2;
t140 = t198 * t125;
t141 = -t126 * pkin(4) - pkin(6) * t171;
t175 = t117 * t125;
t169 = -g(1) + qJDD(2);
t51 = t123 * t169 + t196 * t126;
t94 = t144 * qJD(1);
t42 = t94 * t166 + t51;
t22 = t141 * qJDD(1) + (-t42 + (pkin(6) * t123 * t126 - pkin(4) * t175) * t132) * t122 + t140;
t116 = t122 ^ 2;
t174 = t117 * t132;
t105 = t116 * t174;
t27 = t198 * t122 + t125 * t42;
t89 = t141 * qJD(1);
t23 = -pkin(4) * t105 - pkin(6) * t153 + t89 * t166 + t27;
t10 = t128 * t22 + t130 * t23;
t9 = t128 * t23 - t130 * t22;
t4 = t128 * t10 - t130 * t9;
t188 = t125 * t4;
t172 = t122 * t132;
t149 = t172 * t175;
t83 = -t149 + t162;
t185 = t122 * t83;
t84 = -t149 - t162;
t184 = t125 * t84;
t113 = t126 * t169;
t161 = qJDD(4) - t113;
t168 = t189 + t94;
t29 = -pkin(6) * t105 + (pkin(4) * qJDD(1) * t122 + t64 + (t125 * t89 + t168) * qJD(1)) * t123 + t161;
t183 = t128 * t29;
t54 = t107 - t187;
t182 = t128 * t54;
t181 = t130 * t29;
t180 = t130 * t54;
t178 = t108 * t128;
t177 = t108 * t130;
t118 = t125 ^ 2;
t176 = t117 * t118;
t173 = t119 * t132;
t170 = t126 * t132;
t164 = qJDD(1) * t125;
t160 = t126 * t187;
t158 = t118 * t174;
t157 = t122 * t170;
t155 = t125 * t170;
t5 = t130 * t10 + t128 * t9;
t50 = t196 * t123 - t113;
t28 = t123 * t50 + t126 * t51;
t148 = t122 * t155;
t26 = t122 * t42 - t140;
t14 = t122 * t27 - t125 * t26;
t143 = t117 * t148;
t137 = t144 - t154;
t136 = qJDD(1) * t138;
t135 = t142 * t163;
t115 = t119 * qJDD(1);
t114 = t117 * qJDD(1);
t98 = t123 * t155;
t97 = t115 + t114;
t92 = (-t119 - t176) * t132;
t90 = -t105 - t173;
t85 = t105 + t158;
t82 = (t157 - t164) * t123;
t81 = (t157 + t164) * t123;
t80 = -t98 + t153;
t79 = t98 + t153;
t71 = -t75 + t106;
t70 = t74 - t106;
t69 = -t75 - t106;
t68 = -t122 * t92 + t125 * t83;
t67 = -t122 * t84 + t125 * t90;
t66 = t125 * t92 + t185;
t65 = t122 * t90 + t184;
t61 = t75 - t74;
t59 = -t78 * qJD(5) - t136;
t58 = -t122 * t82 - t125 * t79;
t57 = -t122 * t79 + t125 * t82;
t52 = -t106 - t74;
t46 = t123 * t81 + t126 * t68;
t45 = t123 * t80 + t126 * t67;
t44 = -t74 - t75;
t41 = (t168 * qJD(1) + t64) * t123 + t161;
t40 = t60 - t179;
t37 = -t136 - t194;
t36 = t135 + t194;
t35 = (qJD(5) - t108) * t78 + t135;
t33 = -t128 * t69 + t180;
t32 = t130 * t69 + t182;
t31 = t130 * t52 - t193;
t30 = t128 * t52 + t192;
t25 = t128 * t40 + t130 * t37;
t24 = t128 * t37 - t130 * t40;
t20 = -t122 * t32 + t125 * t33;
t19 = t122 * t33 + t125 * t32;
t18 = -t122 * t30 + t125 * t31;
t17 = t122 * t31 + t125 * t30;
t16 = t123 * t197 + t126 * t20;
t15 = t122 * t26 + t125 * t27;
t13 = t123 * t35 + t126 * t18;
t12 = -t122 * t24 + t125 * t25;
t11 = t122 * t25 + t125 * t24;
t6 = t126 * t12 + t123 * t44;
t3 = -t122 * t4 + t125 * t5;
t2 = t122 * t5 + t188;
t1 = t123 * t29 + t126 * t3;
t7 = [0, 0, 0, 0, 0, qJDD(1), t147, t146, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t127 * qJDD(1) - t124 * t132) + t151, pkin(1) * (-t124 * qJDD(1) - t127 * t132) - t186, 0, pkin(1) * (t124 * t186 + t127 * t151), t114, 0.2e1 * t123 * t162, 0, t115, 0, 0, -t195 * t126, t195 * t123, pkin(2) * t99 + qJ(3) * t97 + pkin(1) * (t124 * t97 + t127 * t99) + t28, -pkin(2) * t63 + qJ(3) * t28 + pkin(1) * (t124 * t28 - t127 * t63), -t143 + (qJDD(1) * t118 + t148) * t117, t123 * (-t122 * t81 - t125 * t80) + t126 * (t105 - t158), t123 * (t184 - (t119 - t176) * t172) + t126 * t82, t143 + (qJDD(1) * t116 - t148) * t117, t123 * (t125 * (t105 - t173) + t185) + t126 * t79, t115, t123 * (-qJ(4) * t65 + t122 * t41) + t126 * (-pkin(3) * t65 + t26) - pkin(2) * t65 + qJ(3) * t45 + pkin(1) * (t124 * t45 - t127 * t65), t123 * (-qJ(4) * t66 + t125 * t41) + t126 * (-pkin(3) * t66 + t27) - pkin(2) * t66 + qJ(3) * t46 + pkin(1) * (t124 * t46 - t127 * t66), -t123 * t14 + t152 * (-t123 * t85 + t126 * t58) + t137 * t57, t152 * (t123 * t41 + t126 * t15) + t137 * t14, t123 * (t125 * (t130 * t60 + t78 * t178) - t122 * (t128 * t60 - t78 * t177)) - t160, t123 * (t125 * (-t128 * t197 - t130 * t35) - t122 * (-t128 * t35 + t130 * t197)) - t126 * t61, t123 * (t125 * (-t128 * t71 + t192) - t122 * (t130 * t71 + t193)) - t126 * t40, t123 * (t125 * (-t128 * t59 - t76 * t177) - t122 * (t130 * t59 - t76 * t178)) + t160, t123 * (t125 * (t130 * t70 + t182) - t122 * (t128 * t70 - t180)) + t126 * t36, t126 * t107 + t123 * (t125 * (-t128 * t78 + t130 * t76) - t122 * (t128 * t76 + t130 * t78)) * t108, t123 * (t125 * (-pkin(6) * t30 + t183) - t122 * (-pkin(4) * t35 + pkin(6) * t31 - t181) - qJ(4) * t17) + t126 * (-pkin(3) * t17 - pkin(4) * t30 + t9) - pkin(2) * t17 + qJ(3) * t13 + pkin(1) * (t124 * t13 - t127 * t17), t123 * (t125 * (-pkin(6) * t32 + t181) - t122 * (-pkin(4) * t197 + pkin(6) * t33 + t183) - qJ(4) * t19) + t126 * (-pkin(3) * t19 - pkin(4) * t32 + t10) - pkin(2) * t19 + qJ(3) * t16 + pkin(1) * (t124 * t16 - t127 * t19), t123 * (t125 * (-pkin(6) * t24 - t4) - t122 * (-pkin(4) * t44 + pkin(6) * t25 + t5) - qJ(4) * t11) + t126 * (-pkin(3) * t11 - pkin(4) * t24) - pkin(2) * t11 + qJ(3) * t6 + pkin(1) * (-t127 * t11 + t124 * t6), t123 * (-pkin(6) * t188 - t122 * (-pkin(4) * t29 + pkin(6) * t5) - qJ(4) * t2) + t126 * (-pkin(3) * t2 - pkin(4) * t4) - pkin(2) * t2 + qJ(3) * t1 + pkin(1) * (t124 * t1 - t127 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123 * t51 - t126 * t50, 0, 0, 0, 0, 0, 0, t123 * t67 - t126 * t80, t123 * t68 - t126 * t81, t123 * t58 + t126 * t85, t123 * t15 - t126 * t41, 0, 0, 0, 0, 0, 0, t123 * t18 - t126 * t35, t123 * t20 - t126 * t197, t123 * t12 - t126 * t44, t123 * t3 - t126 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, t163, -t99, t63, 0, 0, 0, 0, 0, 0, t65, t66, t57, t14, 0, 0, 0, 0, 0, 0, t17, t19, t11, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t81, -t85, t41, 0, 0, 0, 0, 0, 0, t35, t197, t44, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, t61, t40, -t187, -t36, -t107, -t9, -t10, 0, 0;];
tauJ_reg = t7;
