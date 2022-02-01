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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
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
% StartTime: 2022-01-20 09:12:45
% EndTime: 2022-01-20 09:12:51
% DurationCPUTime: 2.19s
% Computational Cost: add. (4582->239), mult. (11021->365), div. (0->0), fcn. (7063->10), ass. (0->158)
t124 = sin(pkin(8));
t118 = t124 ^ 2;
t127 = cos(pkin(8));
t120 = t127 ^ 2;
t133 = qJD(1) ^ 2;
t99 = (t118 + t120) * t133;
t145 = -t127 * pkin(3) - t124 * qJ(4);
t168 = qJD(1) * t124;
t125 = sin(pkin(7));
t128 = cos(pkin(7));
t130 = sin(qJ(1));
t132 = cos(qJ(1));
t152 = t130 * g(1) - t132 * g(2);
t95 = qJDD(1) * pkin(1) + t152;
t147 = t132 * g(1) + t130 * g(2);
t96 = -t133 * pkin(1) - t147;
t151 = -t125 * t96 + t128 * t95;
t63 = -qJDD(1) * pkin(2) - t133 * qJ(3) + qJDD(3) - t151;
t199 = -0.2e1 * qJD(4) * t168 + qJDD(1) * t145 + t63;
t167 = t127 * qJD(1);
t108 = -qJD(5) + t167;
t123 = sin(pkin(9));
t126 = cos(pkin(9));
t129 = sin(qJ(5));
t131 = cos(qJ(5));
t143 = t123 * t131 + t126 * t129;
t139 = t124 * t143;
t76 = qJD(1) * t139;
t180 = t76 * t108;
t164 = t124 * qJDD(1);
t154 = t123 * t164;
t172 = t124 * t126;
t157 = t131 * t172;
t60 = -t76 * qJD(5) + qJDD(1) * t157 - t129 * t154;
t198 = t60 + t180;
t190 = 2 * qJD(3);
t187 = t125 * t95 + t128 * t96;
t64 = -t133 * pkin(2) + qJDD(1) * qJ(3) + t187;
t197 = qJD(1) * t190 + t64;
t153 = pkin(1) * t125 + qJ(3);
t155 = pkin(1) * t128 + pkin(2);
t196 = -t155 * qJDD(1) + t153 * t99 + t63;
t78 = -t129 * t123 * t168 + qJD(1) * t157;
t195 = (qJD(5) + t108) * t78;
t163 = t127 * qJDD(1);
t107 = -qJDD(5) + t163;
t188 = t78 * t76;
t140 = -t107 - t188;
t194 = t129 * t140;
t193 = t131 * t140;
t74 = t76 ^ 2;
t75 = t78 ^ 2;
t106 = t108 ^ 2;
t141 = t199 * t126;
t142 = -t127 * pkin(4) - pkin(6) * t172;
t176 = t118 * t126;
t170 = -g(3) + qJDD(2);
t51 = t124 * t170 + t197 * t127;
t94 = t145 * qJD(1);
t42 = t94 * t167 + t51;
t22 = t142 * qJDD(1) + (-t42 + (pkin(6) * t124 * t127 - pkin(4) * t176) * t133) * t123 + t141;
t117 = t123 ^ 2;
t175 = t118 * t133;
t105 = t117 * t175;
t27 = t199 * t123 + t126 * t42;
t89 = t142 * qJD(1);
t23 = -pkin(4) * t105 - pkin(6) * t154 + t89 * t167 + t27;
t10 = t129 * t22 + t131 * t23;
t9 = t129 * t23 - t131 * t22;
t4 = t129 * t10 - t131 * t9;
t189 = t126 * t4;
t173 = t123 * t133;
t149 = t173 * t176;
t83 = -t149 + t163;
t186 = t123 * t83;
t84 = -t149 - t163;
t185 = t126 * t84;
t113 = t127 * t170;
t162 = qJDD(4) - t113;
t169 = t190 + t94;
t29 = -pkin(6) * t105 + (pkin(4) * qJDD(1) * t123 + t64 + (t126 * t89 + t169) * qJD(1)) * t124 + t162;
t184 = t129 * t29;
t54 = t107 - t188;
t183 = t129 * t54;
t182 = t131 * t29;
t181 = t131 * t54;
t179 = t108 * t129;
t178 = t108 * t131;
t119 = t126 ^ 2;
t177 = t118 * t119;
t174 = t120 * t133;
t171 = t127 * t133;
t165 = qJDD(1) * t126;
t161 = t127 * t188;
t159 = t119 * t175;
t158 = t123 * t171;
t156 = t126 * t171;
t5 = t131 * t10 + t129 * t9;
t50 = t197 * t124 - t113;
t28 = t124 * t50 + t127 * t51;
t148 = t123 * t156;
t26 = t123 * t42 - t141;
t14 = t123 * t27 - t126 * t26;
t144 = t118 * t148;
t138 = t145 - t155;
t137 = qJDD(1) * t139;
t136 = t143 * t164;
t115 = t120 * qJDD(1);
t114 = t118 * qJDD(1);
t98 = t124 * t156;
t97 = t115 + t114;
t92 = (-t120 - t177) * t133;
t90 = -t105 - t174;
t85 = t105 + t159;
t82 = (t158 - t165) * t124;
t81 = (t158 + t165) * t124;
t80 = -t98 + t154;
t79 = t98 + t154;
t71 = -t75 + t106;
t70 = t74 - t106;
t69 = -t75 - t106;
t68 = -t123 * t92 + t126 * t83;
t67 = -t123 * t84 + t126 * t90;
t66 = t126 * t92 + t186;
t65 = t123 * t90 + t185;
t61 = t75 - t74;
t59 = -t78 * qJD(5) - t137;
t58 = -t123 * t82 - t126 * t79;
t57 = -t123 * t79 + t126 * t82;
t52 = -t106 - t74;
t46 = t124 * t81 + t127 * t68;
t45 = t124 * t80 + t127 * t67;
t44 = -t74 - t75;
t41 = (t169 * qJD(1) + t64) * t124 + t162;
t40 = t60 - t180;
t37 = -t137 - t195;
t36 = t136 + t195;
t35 = (qJD(5) - t108) * t78 + t136;
t33 = -t129 * t69 + t181;
t32 = t131 * t69 + t183;
t31 = t131 * t52 - t194;
t30 = t129 * t52 + t193;
t25 = t129 * t40 + t131 * t37;
t24 = t129 * t37 - t131 * t40;
t20 = -t123 * t32 + t126 * t33;
t19 = t123 * t33 + t126 * t32;
t18 = -t123 * t30 + t126 * t31;
t17 = t123 * t31 + t126 * t30;
t16 = t124 * t198 + t127 * t20;
t15 = t123 * t26 + t126 * t27;
t13 = t124 * t35 + t127 * t18;
t12 = -t123 * t24 + t126 * t25;
t11 = t123 * t25 + t126 * t24;
t6 = t127 * t12 + t124 * t44;
t3 = -t123 * t4 + t126 * t5;
t2 = t123 * t5 + t189;
t1 = t124 * t29 + t127 * t3;
t7 = [0, 0, 0, 0, 0, qJDD(1), t152, t147, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t128 * qJDD(1) - t125 * t133) + t151, pkin(1) * (-t125 * qJDD(1) - t128 * t133) - t187, 0, pkin(1) * (t125 * t187 + t128 * t151), t114, 0.2e1 * t124 * t163, 0, t115, 0, 0, -t196 * t127, t196 * t124, pkin(2) * t99 + qJ(3) * t97 + pkin(1) * (t125 * t97 + t128 * t99) + t28, -pkin(2) * t63 + qJ(3) * t28 + pkin(1) * (t125 * t28 - t128 * t63), -t144 + (qJDD(1) * t119 + t148) * t118, t124 * (-t123 * t81 - t126 * t80) + t127 * (t105 - t159), t124 * (t185 - (t120 - t177) * t173) + t127 * t82, t144 + (qJDD(1) * t117 - t148) * t118, t124 * (t126 * (t105 - t174) + t186) + t127 * t79, t115, t124 * (-qJ(4) * t65 + t123 * t41) + t127 * (-pkin(3) * t65 + t26) - pkin(2) * t65 + qJ(3) * t45 + pkin(1) * (t125 * t45 - t128 * t65), t124 * (-qJ(4) * t66 + t126 * t41) + t127 * (-pkin(3) * t66 + t27) - pkin(2) * t66 + qJ(3) * t46 + pkin(1) * (t125 * t46 - t128 * t66), -t124 * t14 + t153 * (-t124 * t85 + t127 * t58) + t138 * t57, t153 * (t124 * t41 + t127 * t15) + t138 * t14, t124 * (t126 * (t131 * t60 + t78 * t179) - t123 * (t129 * t60 - t78 * t178)) - t161, t124 * (t126 * (-t129 * t198 - t131 * t35) - t123 * (-t129 * t35 + t131 * t198)) - t127 * t61, t124 * (t126 * (-t129 * t71 + t193) - t123 * (t131 * t71 + t194)) - t127 * t40, t124 * (t126 * (-t129 * t59 - t76 * t178) - t123 * (t131 * t59 - t76 * t179)) + t161, t124 * (t126 * (t131 * t70 + t183) - t123 * (t129 * t70 - t181)) + t127 * t36, t127 * t107 + t124 * (t126 * (-t129 * t78 + t131 * t76) - t123 * (t129 * t76 + t131 * t78)) * t108, t124 * (t126 * (-pkin(6) * t30 + t184) - t123 * (-pkin(4) * t35 + pkin(6) * t31 - t182) - qJ(4) * t17) + t127 * (-pkin(3) * t17 - pkin(4) * t30 + t9) - pkin(2) * t17 + qJ(3) * t13 + pkin(1) * (t125 * t13 - t128 * t17), t124 * (t126 * (-pkin(6) * t32 + t182) - t123 * (-pkin(4) * t198 + pkin(6) * t33 + t184) - qJ(4) * t19) + t127 * (-pkin(3) * t19 - pkin(4) * t32 + t10) - pkin(2) * t19 + qJ(3) * t16 + pkin(1) * (t125 * t16 - t128 * t19), t124 * (t126 * (-pkin(6) * t24 - t4) - t123 * (-pkin(4) * t44 + pkin(6) * t25 + t5) - qJ(4) * t11) + t127 * (-pkin(3) * t11 - pkin(4) * t24) - pkin(2) * t11 + qJ(3) * t6 + pkin(1) * (-t128 * t11 + t125 * t6), t124 * (-pkin(6) * t189 - t123 * (-pkin(4) * t29 + pkin(6) * t5) - qJ(4) * t2) + t127 * (-pkin(3) * t2 - pkin(4) * t4) - pkin(2) * t2 + qJ(3) * t1 + pkin(1) * (t125 * t1 - t128 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124 * t51 - t127 * t50, 0, 0, 0, 0, 0, 0, t124 * t67 - t127 * t80, t124 * t68 - t127 * t81, t124 * t58 + t127 * t85, t124 * t15 - t127 * t41, 0, 0, 0, 0, 0, 0, t124 * t18 - t127 * t35, t124 * t20 - t127 * t198, t124 * t12 - t127 * t44, t124 * t3 - t127 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, t164, -t99, t63, 0, 0, 0, 0, 0, 0, t65, t66, t57, t14, 0, 0, 0, 0, 0, 0, t17, t19, t11, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t81, -t85, t41, 0, 0, 0, 0, 0, 0, t35, t198, t44, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, t61, t40, -t188, -t36, -t107, -t9, -t10, 0, 0;];
tauJ_reg = t7;
