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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:29:07
% EndTime: 2019-12-05 17:29:15
% DurationCPUTime: 2.00s
% Computational Cost: add. (4582->239), mult. (11021->365), div. (0->0), fcn. (7063->10), ass. (0->158)
t125 = sin(pkin(8));
t119 = t125 ^ 2;
t128 = cos(pkin(8));
t121 = t128 ^ 2;
t134 = qJD(1) ^ 2;
t99 = (t119 + t121) * t134;
t146 = -t128 * pkin(3) - t125 * qJ(4);
t168 = qJD(1) * t125;
t126 = sin(pkin(7));
t129 = cos(pkin(7));
t131 = sin(qJ(1));
t133 = cos(qJ(1));
t169 = t133 * g(2) + t131 * g(3);
t95 = qJDD(1) * pkin(1) + t169;
t148 = -t131 * g(2) + t133 * g(3);
t96 = -t134 * pkin(1) - t148;
t152 = -t126 * t96 + t129 * t95;
t63 = -qJDD(1) * pkin(2) - t134 * qJ(3) + qJDD(3) - t152;
t200 = -0.2e1 * qJD(4) * t168 + qJDD(1) * t146 + t63;
t167 = t128 * qJD(1);
t108 = -qJD(5) + t167;
t124 = sin(pkin(9));
t127 = cos(pkin(9));
t130 = sin(qJ(5));
t132 = cos(qJ(5));
t144 = t124 * t132 + t127 * t130;
t140 = t125 * t144;
t76 = qJD(1) * t140;
t181 = t76 * t108;
t164 = t125 * qJDD(1);
t154 = t124 * t164;
t173 = t125 * t127;
t157 = t132 * t173;
t60 = -t76 * qJD(5) + qJDD(1) * t157 - t130 * t154;
t199 = t60 + t181;
t191 = 2 * qJD(3);
t188 = t126 * t95 + t129 * t96;
t64 = -t134 * pkin(2) + qJDD(1) * qJ(3) + t188;
t198 = qJD(1) * t191 + t64;
t153 = pkin(1) * t126 + qJ(3);
t155 = pkin(1) * t129 + pkin(2);
t197 = -t155 * qJDD(1) + t153 * t99 + t63;
t78 = -t130 * t124 * t168 + qJD(1) * t157;
t196 = (qJD(5) + t108) * t78;
t163 = t128 * qJDD(1);
t107 = -qJDD(5) + t163;
t189 = t78 * t76;
t141 = -t107 - t189;
t195 = t130 * t141;
t194 = t132 * t141;
t74 = t76 ^ 2;
t75 = t78 ^ 2;
t106 = t108 ^ 2;
t142 = t200 * t127;
t143 = -t128 * pkin(4) - pkin(6) * t173;
t177 = t119 * t127;
t171 = -g(1) + qJDD(2);
t51 = t125 * t171 + t198 * t128;
t94 = t146 * qJD(1);
t42 = t94 * t167 + t51;
t22 = t143 * qJDD(1) + (-t42 + (pkin(6) * t125 * t128 - pkin(4) * t177) * t134) * t124 + t142;
t118 = t124 ^ 2;
t176 = t119 * t134;
t105 = t118 * t176;
t27 = t200 * t124 + t127 * t42;
t89 = t143 * qJD(1);
t23 = -pkin(4) * t105 - pkin(6) * t154 + t89 * t167 + t27;
t10 = t130 * t22 + t132 * t23;
t9 = t130 * t23 - t132 * t22;
t4 = t130 * t10 - t132 * t9;
t190 = t127 * t4;
t174 = t124 * t134;
t150 = t174 * t177;
t83 = -t150 + t163;
t187 = t124 * t83;
t84 = -t150 - t163;
t186 = t127 * t84;
t113 = t128 * t171;
t162 = qJDD(4) - t113;
t170 = t191 + t94;
t29 = -pkin(6) * t105 + (pkin(4) * qJDD(1) * t124 + t64 + (t127 * t89 + t170) * qJD(1)) * t125 + t162;
t185 = t130 * t29;
t54 = t107 - t189;
t184 = t130 * t54;
t183 = t132 * t29;
t182 = t132 * t54;
t180 = t108 * t130;
t179 = t108 * t132;
t120 = t127 ^ 2;
t178 = t119 * t120;
t175 = t121 * t134;
t172 = t128 * t134;
t165 = qJDD(1) * t127;
t161 = t128 * t189;
t159 = t120 * t176;
t158 = t124 * t172;
t156 = t127 * t172;
t5 = t132 * t10 + t130 * t9;
t50 = t198 * t125 - t113;
t28 = t125 * t50 + t128 * t51;
t149 = t124 * t156;
t26 = t124 * t42 - t142;
t14 = t124 * t27 - t127 * t26;
t145 = t119 * t149;
t139 = t146 - t155;
t138 = qJDD(1) * t140;
t137 = t144 * t164;
t115 = t121 * qJDD(1);
t114 = t119 * qJDD(1);
t98 = t125 * t156;
t97 = t115 + t114;
t92 = (-t121 - t178) * t134;
t90 = -t105 - t175;
t85 = t105 + t159;
t82 = (t158 - t165) * t125;
t81 = (t158 + t165) * t125;
t80 = -t98 + t154;
t79 = t98 + t154;
t71 = -t75 + t106;
t70 = t74 - t106;
t69 = -t75 - t106;
t68 = -t124 * t92 + t127 * t83;
t67 = -t124 * t84 + t127 * t90;
t66 = t127 * t92 + t187;
t65 = t124 * t90 + t186;
t61 = t75 - t74;
t59 = -t78 * qJD(5) - t138;
t58 = -t124 * t82 - t127 * t79;
t57 = -t124 * t79 + t127 * t82;
t52 = -t106 - t74;
t46 = t125 * t81 + t128 * t68;
t45 = t125 * t80 + t128 * t67;
t44 = -t74 - t75;
t41 = (t170 * qJD(1) + t64) * t125 + t162;
t40 = t60 - t181;
t37 = -t138 - t196;
t36 = t137 + t196;
t35 = (qJD(5) - t108) * t78 + t137;
t33 = -t130 * t69 + t182;
t32 = t132 * t69 + t184;
t31 = t132 * t52 - t195;
t30 = t130 * t52 + t194;
t25 = t130 * t40 + t132 * t37;
t24 = t130 * t37 - t132 * t40;
t20 = -t124 * t32 + t127 * t33;
t19 = t124 * t33 + t127 * t32;
t18 = -t124 * t30 + t127 * t31;
t17 = t124 * t31 + t127 * t30;
t16 = t125 * t199 + t128 * t20;
t15 = t124 * t26 + t127 * t27;
t13 = t125 * t35 + t128 * t18;
t12 = -t124 * t24 + t127 * t25;
t11 = t124 * t25 + t127 * t24;
t6 = t128 * t12 + t125 * t44;
t3 = -t124 * t4 + t127 * t5;
t2 = t124 * t5 + t190;
t1 = t125 * t29 + t128 * t3;
t7 = [0, 0, 0, 0, 0, qJDD(1), t169, t148, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t129 * qJDD(1) - t126 * t134) + t152, pkin(1) * (-t126 * qJDD(1) - t129 * t134) - t188, 0, pkin(1) * (t126 * t188 + t129 * t152), t114, 0.2e1 * t125 * t163, 0, t115, 0, 0, -t197 * t128, t197 * t125, pkin(2) * t99 + qJ(3) * t97 + pkin(1) * (t126 * t97 + t129 * t99) + t28, -pkin(2) * t63 + qJ(3) * t28 + pkin(1) * (t126 * t28 - t129 * t63), -t145 + (qJDD(1) * t120 + t149) * t119, t125 * (-t124 * t81 - t127 * t80) + t128 * (t105 - t159), t125 * (t186 - (t121 - t178) * t174) + t128 * t82, t145 + (qJDD(1) * t118 - t149) * t119, t125 * (t127 * (t105 - t175) + t187) + t128 * t79, t115, t125 * (-qJ(4) * t65 + t124 * t41) + t128 * (-pkin(3) * t65 + t26) - pkin(2) * t65 + qJ(3) * t45 + pkin(1) * (t126 * t45 - t129 * t65), t125 * (-qJ(4) * t66 + t127 * t41) + t128 * (-pkin(3) * t66 + t27) - pkin(2) * t66 + qJ(3) * t46 + pkin(1) * (t126 * t46 - t129 * t66), -t125 * t14 + t153 * (-t125 * t85 + t128 * t58) + t139 * t57, t153 * (t125 * t41 + t128 * t15) + t139 * t14, t125 * (t127 * (t132 * t60 + t78 * t180) - t124 * (t130 * t60 - t78 * t179)) - t161, t125 * (t127 * (-t130 * t199 - t132 * t35) - t124 * (-t130 * t35 + t132 * t199)) - t128 * t61, t125 * (t127 * (-t130 * t71 + t194) - t124 * (t132 * t71 + t195)) - t128 * t40, t125 * (t127 * (-t130 * t59 - t76 * t179) - t124 * (t132 * t59 - t76 * t180)) + t161, t125 * (t127 * (t132 * t70 + t184) - t124 * (t130 * t70 - t182)) + t128 * t36, t128 * t107 + t125 * (t127 * (-t130 * t78 + t132 * t76) - t124 * (t130 * t76 + t132 * t78)) * t108, t125 * (t127 * (-pkin(6) * t30 + t185) - t124 * (-pkin(4) * t35 + pkin(6) * t31 - t183) - qJ(4) * t17) + t128 * (-pkin(3) * t17 - pkin(4) * t30 + t9) - pkin(2) * t17 + qJ(3) * t13 + pkin(1) * (t126 * t13 - t129 * t17), t125 * (t127 * (-pkin(6) * t32 + t183) - t124 * (-pkin(4) * t199 + pkin(6) * t33 + t185) - qJ(4) * t19) + t128 * (-pkin(3) * t19 - pkin(4) * t32 + t10) - pkin(2) * t19 + qJ(3) * t16 + pkin(1) * (t126 * t16 - t129 * t19), t125 * (t127 * (-pkin(6) * t24 - t4) - t124 * (-pkin(4) * t44 + pkin(6) * t25 + t5) - qJ(4) * t11) + t128 * (-pkin(3) * t11 - pkin(4) * t24) - pkin(2) * t11 + qJ(3) * t6 + pkin(1) * (-t129 * t11 + t126 * t6), t125 * (-pkin(6) * t190 - t124 * (-pkin(4) * t29 + pkin(6) * t5) - qJ(4) * t2) + t128 * (-pkin(3) * t2 - pkin(4) * t4) - pkin(2) * t2 + qJ(3) * t1 + pkin(1) * (t126 * t1 - t129 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125 * t51 - t128 * t50, 0, 0, 0, 0, 0, 0, t125 * t67 - t128 * t80, t125 * t68 - t128 * t81, t125 * t58 + t128 * t85, t125 * t15 - t128 * t41, 0, 0, 0, 0, 0, 0, t125 * t18 - t128 * t35, t125 * t20 - t128 * t199, t125 * t12 - t128 * t44, t125 * t3 - t128 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, t164, -t99, t63, 0, 0, 0, 0, 0, 0, t65, t66, t57, t14, 0, 0, 0, 0, 0, 0, t17, t19, t11, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t81, -t85, t41, 0, 0, 0, 0, 0, 0, t35, t199, t44, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t61, t40, -t189, -t36, -t107, -t9, -t10, 0, 0;];
tauJ_reg = t7;
