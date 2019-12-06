% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRR1
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:07
% EndTime: 2019-12-05 15:43:14
% DurationCPUTime: 2.05s
% Computational Cost: add. (5858->246), mult. (13688->367), div. (0->0), fcn. (10232->10), ass. (0->157)
t166 = sin(pkin(8));
t167 = cos(pkin(8));
t114 = -t167 * g(1) - t166 * g(2);
t135 = sin(qJ(2));
t147 = t166 * g(1) - t167 * g(2);
t180 = cos(qJ(2));
t144 = -t180 * t114 - t135 * t147;
t181 = qJD(2) ^ 2;
t190 = -t181 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) - t144;
t134 = sin(qJ(4));
t131 = sin(pkin(9));
t132 = cos(pkin(9));
t137 = cos(qJ(4));
t109 = (-t131 * t134 + t132 * t137) * qJD(2);
t148 = t131 * t137 + t132 * t134;
t111 = t148 * qJD(2);
t165 = t111 * t109;
t187 = qJDD(4) + t165;
t192 = t134 * t187;
t191 = t137 * t187;
t133 = sin(qJ(5));
t126 = qJDD(4) + qJDD(5);
t136 = cos(qJ(5));
t82 = -t136 * t109 + t133 * t111;
t84 = t133 * t109 + t136 * t111;
t61 = t84 * t82;
t185 = -t61 + t126;
t189 = t133 * t185;
t188 = t136 * t185;
t128 = t181 * qJ(3);
t130 = qJDD(2) * pkin(2);
t139 = t131 ^ 2;
t141 = t132 ^ 2;
t184 = t139 + t141;
t150 = -t135 * t114 + t180 * t147;
t92 = qJDD(3) - t128 - t130 - t150;
t186 = t184 * t128 - t130 + t92;
t129 = -g(3) + qJDD(1);
t119 = t132 * t129;
t183 = t119 + (pkin(3) * t181 * t132 - pkin(6) * qJDD(2) - t190) * t131;
t161 = t109 * qJD(4);
t155 = t132 * qJDD(2);
t159 = t141 * t181;
t74 = t131 * t129 + t190 * t132;
t67 = -pkin(3) * t159 + pkin(6) * t155 + t74;
t45 = t134 * t67 - t137 * t183;
t108 = t148 * qJDD(2);
t96 = t108 + t161;
t182 = -t45 + (-t96 + t161) * pkin(7);
t80 = t82 ^ 2;
t81 = t84 ^ 2;
t106 = t109 ^ 2;
t107 = t111 ^ 2;
t127 = qJD(4) + qJD(5);
t125 = t127 ^ 2;
t143 = pkin(4) * t187 + t182;
t46 = t183 * t134 + t137 * t67;
t160 = t111 * qJD(4);
t156 = t131 * qJDD(2);
t72 = -t134 * t156 + t137 * t155;
t94 = t72 - t160;
t99 = qJD(4) * pkin(4) - t111 * pkin(7);
t28 = -t106 * pkin(4) + t94 * pkin(7) - qJD(4) * t99 + t46;
t13 = t133 * t28 - t136 * t143;
t172 = t136 * t28;
t14 = t133 * t143 + t172;
t6 = -t136 * t13 + t133 * t14;
t179 = t134 * t6;
t178 = t137 * t6;
t25 = t134 * t46 - t137 * t45;
t177 = t131 * t25;
t71 = -pkin(3) * t155 + t92 + (-t139 * t181 - t159) * pkin(6);
t41 = -t94 * pkin(4) - t106 * pkin(7) + t111 * t99 + t71;
t176 = t133 * t41;
t58 = t61 + t126;
t175 = t133 * t58;
t174 = t134 * t71;
t90 = qJDD(4) - t165;
t173 = t134 * t90;
t171 = t136 * t41;
t170 = t136 * t58;
t169 = t137 * t71;
t168 = t137 * t90;
t164 = t127 * t133;
t163 = t127 * t136;
t158 = qJD(5) + t127;
t7 = t133 * t13 + t136 * t14;
t73 = t190 * t131 - t119;
t153 = t131 * t73 + t132 * t74;
t152 = t133 * t96 - t136 * t94;
t26 = t134 * t45 + t137 * t46;
t149 = t133 * t94 + t136 * t96;
t146 = (-qJD(5) + t127) * t84 - t152;
t53 = -t82 * qJD(5) + t149;
t138 = qJD(4) ^ 2;
t124 = t141 * qJDD(2);
t123 = t139 * qJDD(2);
t113 = t184 * t181;
t102 = -t107 - t138;
t101 = -t107 + t138;
t100 = t106 - t138;
t95 = t108 + 0.2e1 * t161;
t93 = -t72 + 0.2e1 * t160;
t88 = -t138 - t106;
t79 = t127 * t82;
t78 = -t81 + t125;
t77 = t80 - t125;
t76 = -t81 - t125;
t75 = -t106 - t107;
t69 = -t134 * t102 - t168;
t68 = t137 * t102 - t173;
t65 = t134 * t108 + t137 * t72;
t64 = -t137 * t108 + t134 * t72;
t63 = t137 * t88 - t192;
t62 = t134 * t88 + t191;
t60 = t81 - t80;
t56 = -t125 - t80;
t55 = (t133 * t84 - t136 * t82) * t127;
t54 = (-t133 * t82 - t136 * t84) * t127;
t52 = -t84 * qJD(5) - t152;
t51 = -t80 - t81;
t50 = t136 * t77 - t175;
t49 = -t133 * t78 + t188;
t48 = t133 * t77 + t170;
t47 = t136 * t78 + t189;
t44 = -t133 * t76 - t170;
t43 = t136 * t76 - t175;
t40 = t53 + t79;
t39 = t53 - t79;
t38 = -t158 * t82 + t149;
t35 = t158 * t84 + t152;
t34 = t136 * t53 - t84 * t164;
t33 = t133 * t53 + t84 * t163;
t32 = -t133 * t52 + t82 * t163;
t31 = t136 * t52 + t82 * t164;
t30 = t136 * t56 - t189;
t29 = t133 * t56 + t188;
t24 = -t134 * t43 + t137 * t44;
t23 = t134 * t44 + t137 * t43;
t22 = -pkin(7) * t43 + t171;
t21 = -pkin(7) * t29 + t176;
t20 = t133 * t40 + t136 * t146;
t19 = -t133 * t39 - t136 * t35;
t18 = t133 * t146 - t136 * t40;
t17 = -t133 * t35 + t136 * t39;
t16 = -t134 * t29 + t137 * t30;
t15 = t134 * t30 + t137 * t29;
t11 = -pkin(4) * t38 + pkin(7) * t44 + t176;
t10 = -pkin(4) * t35 + pkin(7) * t30 - t171;
t9 = -t134 * t18 + t137 * t20;
t8 = t134 * t20 + t137 * t18;
t5 = -pkin(4) * t41 + pkin(7) * t7;
t4 = -pkin(7) * t18 - t6;
t3 = -pkin(4) * t51 + pkin(7) * t20 + t7;
t2 = t137 * t7 - t179;
t1 = t134 * t7 + t178;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t129, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131 * t74 - t132 * t73, 0, 0, 0, 0, 0, 0, t131 * t63 + t132 * t62, t131 * t69 + t132 * t68, t131 * t65 + t132 * t64, t131 * t26 + t132 * t25, 0, 0, 0, 0, 0, 0, t131 * t16 + t132 * t15, t131 * t24 + t132 * t23, t131 * t9 + t132 * t8, t132 * t1 + t131 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t150, t144, 0, 0, t123, 0.2e1 * t131 * t155, 0, t124, 0, 0, -t186 * t132, t186 * t131, pkin(2) * t113 + qJ(3) * (t124 + t123) + t153, -pkin(2) * t92 + qJ(3) * t153, t131 * (-t134 * t160 + t137 * t96) + t132 * (t134 * t96 + t137 * t160), t131 * (-t134 * t95 - t137 * t93) + t132 * (-t134 * t93 + t137 * t95), t131 * (-t134 * t101 + t191) + t132 * (t137 * t101 + t192), t131 * (-t134 * t94 - t137 * t161) + t132 * (-t134 * t161 + t137 * t94), t131 * (t137 * t100 - t173) + t132 * (t134 * t100 + t168), (t131 * (t109 * t137 + t111 * t134) + t132 * (t109 * t134 - t111 * t137)) * qJD(4), t131 * (-pkin(6) * t62 + t174) + t132 * (-pkin(3) * t93 + pkin(6) * t63 - t169) - pkin(2) * t93 + qJ(3) * (-t131 * t62 + t132 * t63), t131 * (-pkin(6) * t68 + t169) + t132 * (-pkin(3) * t95 + pkin(6) * t69 + t174) - pkin(2) * t95 + qJ(3) * (-t131 * t68 + t132 * t69), t131 * (-pkin(6) * t64 - t25) + t132 * (-pkin(3) * t75 + pkin(6) * t65 + t26) - pkin(2) * t75 + qJ(3) * (-t131 * t64 + t132 * t65), -pkin(6) * t177 + t132 * (-pkin(3) * t71 + pkin(6) * t26) - pkin(2) * t71 + qJ(3) * (t132 * t26 - t177), t131 * (-t134 * t33 + t137 * t34) + t132 * (t134 * t34 + t137 * t33), t131 * (-t134 * t17 + t137 * t19) + t132 * (t134 * t19 + t137 * t17), t131 * (-t134 * t47 + t137 * t49) + t132 * (t134 * t49 + t137 * t47), t131 * (-t134 * t31 + t137 * t32) + t132 * (t134 * t32 + t137 * t31), t131 * (-t134 * t48 + t137 * t50) + t132 * (t134 * t50 + t137 * t48), t131 * (-t134 * t54 + t137 * t55) + t132 * (t134 * t55 + t137 * t54), t131 * (-pkin(6) * t15 - t134 * t10 + t137 * t21) + t132 * (-pkin(3) * t35 + pkin(6) * t16 + t137 * t10 + t134 * t21) - pkin(2) * t35 + qJ(3) * (-t131 * t15 + t132 * t16), t131 * (-pkin(6) * t23 - t134 * t11 + t137 * t22) + t132 * (-pkin(3) * t38 + pkin(6) * t24 + t137 * t11 + t134 * t22) - pkin(2) * t38 + qJ(3) * (-t131 * t23 + t132 * t24), t131 * (-pkin(6) * t8 - t134 * t3 + t137 * t4) + t132 * (-pkin(3) * t51 + pkin(6) * t9 + t134 * t4 + t137 * t3) - pkin(2) * t51 + qJ(3) * (-t131 * t8 + t132 * t9), t131 * (-pkin(6) * t1 - pkin(7) * t178 - t134 * t5) + t132 * (-pkin(3) * t41 + pkin(6) * t2 - pkin(7) * t179 + t137 * t5) - pkin(2) * t41 + qJ(3) * (-t131 * t1 + t132 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, t156, -t113, t92, 0, 0, 0, 0, 0, 0, t93, t95, t75, t71, 0, 0, 0, 0, 0, 0, t35, t38, t51, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, t107 - t106, t108, t165, t72, qJDD(4), -t45, -t46, 0, 0, t61, t60, t40, -t61, t146, t126, pkin(4) * t29 - t13, -t172 - t133 * t182 + (-t133 * t187 + t43) * pkin(4), pkin(4) * t18, pkin(4) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t60, t40, -t61, t146, t126, -t13, -t14, 0, 0;];
tauJ_reg = t12;
