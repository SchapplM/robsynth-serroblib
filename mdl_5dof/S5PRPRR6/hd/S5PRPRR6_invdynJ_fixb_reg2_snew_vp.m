% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:55
% EndTime: 2019-12-05 15:58:01
% DurationCPUTime: 2.10s
% Computational Cost: add. (6357->280), mult. (14180->425), div. (0->0), fcn. (10794->12), ass. (0->166)
t143 = sin(qJ(5));
t138 = sin(pkin(10));
t140 = cos(pkin(10));
t144 = sin(qJ(4));
t147 = cos(qJ(4));
t161 = t138 * t147 + t140 * t144;
t115 = t161 * qJD(2);
t146 = cos(qJ(5));
t101 = t143 * qJD(4) + t146 * t115;
t99 = -t146 * qJD(4) + t143 * t115;
t78 = t101 * t99;
t170 = t140 * qJDD(2);
t171 = t138 * qJDD(2);
t162 = t144 * t171 - t147 * t170;
t175 = t115 * qJD(4);
t93 = -t162 - t175;
t86 = qJDD(5) - t93;
t200 = -t78 + t86;
t206 = t143 * t200;
t177 = qJD(2) * t140;
t180 = t138 * t144;
t113 = qJD(2) * t180 - t147 * t177;
t96 = t115 * t113;
t197 = qJDD(4) - t96;
t205 = t144 * t197;
t204 = t146 * t200;
t203 = t147 * t197;
t133 = t138 ^ 2;
t134 = t140 ^ 2;
t196 = qJD(2) ^ 2;
t120 = (t133 + t134) * t196;
t139 = sin(pkin(5));
t141 = cos(pkin(5));
t183 = sin(pkin(9));
t184 = cos(pkin(9));
t159 = t183 * g(1) - t184 * g(2);
t158 = t141 * t159;
t178 = -g(3) + qJDD(1);
t202 = t139 * t178 + t158;
t157 = -t139 * t159 + t141 * t178;
t155 = t140 * t157;
t201 = qJ(3) + pkin(7);
t121 = -t184 * g(1) - t183 * g(2);
t145 = sin(qJ(2));
t148 = cos(qJ(2));
t85 = t148 * t121 + t202 * t145;
t152 = t155 + (-t201 * qJDD(2) + (-(2 * qJD(3)) + (t140 * pkin(3) + pkin(2)) * qJD(2)) * qJD(2) - t85) * t138;
t174 = t134 * t196;
t154 = qJDD(2) * qJ(3) + t85;
t195 = 2 * qJD(3);
t62 = t140 * (-t196 * pkin(2) + t154) + t138 * t157 + t177 * t195;
t55 = -pkin(3) * t174 + pkin(7) * t170 + t62;
t33 = t144 * t152 + t147 * t55;
t109 = qJD(5) + t113;
t112 = t161 * qJDD(2);
t176 = t113 * qJD(4);
t95 = t112 - t176;
t165 = -t146 * qJDD(4) + t143 * t95;
t48 = (qJD(5) - t109) * t101 + t165;
t97 = t99 ^ 2;
t98 = t101 ^ 2;
t108 = t109 ^ 2;
t110 = t113 ^ 2;
t111 = t115 ^ 2;
t32 = t144 * t55 - t147 * t152;
t18 = t144 * t33 - t147 * t32;
t194 = t138 * t18;
t149 = qJD(4) ^ 2;
t87 = t113 * pkin(4) - t115 * pkin(8);
t26 = -qJDD(4) * pkin(4) - t149 * pkin(8) + t115 * t87 + t32;
t193 = t143 * t26;
t58 = t78 + t86;
t192 = t143 * t58;
t137 = qJDD(2) * pkin(2);
t163 = t145 * t121 - t202 * t148;
t80 = -t196 * qJ(3) + qJDD(3) - t137 + t163;
t71 = -pkin(3) * t170 + t80 + (-t133 * t196 - t174) * pkin(7);
t191 = t144 * t71;
t90 = qJDD(4) + t96;
t190 = t144 * t90;
t189 = t146 * t26;
t188 = t146 * t58;
t186 = t147 * t71;
t185 = t147 * t90;
t182 = t109 * t143;
t181 = t109 * t146;
t172 = qJD(5) + t109;
t169 = t144 * t78;
t168 = t147 * t78;
t167 = -pkin(4) * t147 - pkin(3);
t27 = -t149 * pkin(4) + qJDD(4) * pkin(8) - t113 * t87 + t33;
t35 = (-t95 + t176) * pkin(8) + (-t93 + t175) * pkin(4) + t71;
t14 = t143 * t27 - t146 * t35;
t15 = t143 * t35 + t146 * t27;
t6 = t143 * t14 + t146 * t15;
t61 = -t155 + ((-pkin(2) * qJD(2) + t195) * qJD(2) + t154) * t138;
t36 = t138 * t61 + t140 * t62;
t19 = t144 * t32 + t147 * t33;
t166 = -t80 + t137;
t2 = t144 * t6 - t147 * t26;
t3 = t144 * t26 + t147 * t6;
t1 = -t138 * t2 + t140 * t3;
t5 = -t146 * t14 + t143 * t15;
t160 = -t143 * qJDD(4) - t146 * t95;
t69 = -t99 * qJD(5) - t160;
t130 = t134 * qJDD(2);
t129 = t133 * qJDD(2);
t119 = t130 + t129;
t118 = t140 * t120;
t117 = t138 * t120;
t105 = -t111 - t149;
t104 = -t111 + t149;
t103 = t110 - t149;
t94 = t112 - 0.2e1 * t176;
t92 = t162 + 0.2e1 * t175;
t88 = -t149 - t110;
t83 = t109 * t99;
t82 = -t98 + t108;
t81 = t97 - t108;
t77 = -t110 - t111;
t76 = t98 - t97;
t74 = -t98 - t108;
t73 = -t144 * t105 - t185;
t72 = t147 * t105 - t190;
t70 = -t108 - t97;
t68 = -t101 * qJD(5) - t165;
t67 = t97 + t98;
t66 = t144 * t112 - t147 * t162;
t65 = -t147 * t112 - t144 * t162;
t64 = t147 * t88 - t205;
t63 = t144 * t88 + t203;
t56 = (t101 * t143 - t146 * t99) * t109;
t53 = t172 * t99 + t160;
t52 = t69 + t83;
t51 = t69 - t83;
t50 = -t172 * t101 - t165;
t47 = -t101 * t182 + t146 * t69;
t46 = -t143 * t68 + t99 * t181;
t45 = -t138 * t72 + t140 * t73;
t44 = t146 * t81 - t192;
t43 = -t143 * t82 + t204;
t42 = -t143 * t74 - t188;
t41 = t146 * t74 - t192;
t40 = -t138 * t65 + t140 * t66;
t39 = t146 * t70 - t206;
t38 = t143 * t70 + t204;
t37 = -t138 * t63 + t140 * t64;
t30 = -t143 * t51 + t146 * t50;
t29 = t143 * t52 - t146 * t48;
t28 = -t143 * t48 - t146 * t52;
t25 = -t144 * t53 + t147 * t42;
t24 = t144 * t42 + t147 * t53;
t23 = -t144 * t50 + t147 * t39;
t22 = t144 * t39 + t147 * t50;
t21 = -t144 * t67 + t147 * t29;
t20 = t144 * t29 + t147 * t67;
t17 = -pkin(8) * t41 + t189;
t16 = -pkin(8) * t38 + t193;
t12 = -t138 * t24 + t140 * t25;
t11 = -t138 * t22 + t140 * t23;
t10 = -pkin(4) * t41 + t15;
t9 = -pkin(4) * t38 + t14;
t8 = -t138 * t20 + t140 * t21;
t7 = t140 * t19 - t194;
t4 = -pkin(8) * t28 - t5;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t178, 0, 0, 0, 0, 0, 0, (qJDD(2) * t148 - t196 * t145) * t139, (-qJDD(2) * t145 - t196 * t148) * t139, 0, t141 ^ 2 * t178 + (t145 * t85 - t148 * t163 - t158) * t139, 0, 0, 0, 0, 0, 0, (-t118 * t145 + t148 * t170) * t139, (t117 * t145 - t148 * t171) * t139, (t119 * t145 + t120 * t148) * t139, t141 * (t138 * t62 - t140 * t61) + (t145 * t36 - t148 * t80) * t139, 0, 0, 0, 0, 0, 0, t141 * (t138 * t64 + t140 * t63) + (t145 * t37 - t148 * t92) * t139, t141 * (t138 * t73 + t140 * t72) + (t145 * t45 - t148 * t94) * t139, t141 * (t138 * t66 + t140 * t65) + (t145 * t40 - t148 * t77) * t139, t141 * (t138 * t19 + t140 * t18) + (t145 * t7 - t148 * t71) * t139, 0, 0, 0, 0, 0, 0, t141 * (t138 * t23 + t140 * t22) + (t11 * t145 - t148 * t38) * t139, t141 * (t138 * t25 + t140 * t24) + (t12 * t145 - t148 * t41) * t139, t141 * (t138 * t21 + t140 * t20) + (t145 * t8 - t148 * t28) * t139, t141 * (t138 * t3 + t140 * t2) + (t1 * t145 - t148 * t5) * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t163, -t85, 0, 0, t129, 0.2e1 * t138 * t170, 0, t130, 0, 0, -qJ(3) * t118 + t166 * t140, qJ(3) * t117 - t166 * t138, pkin(2) * t120 + qJ(3) * t119 + t36, -pkin(2) * t80 + qJ(3) * t36, t138 * (-t144 * t175 + t147 * t95) + t140 * (t144 * t95 + t147 * t175), t138 * (-t144 * t94 - t147 * t92) + t140 * (-t144 * t92 + t147 * t94), t138 * (-t104 * t144 + t203) + t140 * (t104 * t147 + t205), t138 * (-t144 * t93 + t147 * t176) + t140 * (t144 * t176 + t147 * t93), t138 * (t103 * t147 - t190) + t140 * (t103 * t144 + t185), (t138 * (-t113 * t147 + t115 * t144) + t140 * (-t113 * t144 - t115 * t147)) * qJD(4), t138 * (-pkin(7) * t63 + t191) + t140 * (-pkin(3) * t92 + pkin(7) * t64 - t186) - pkin(2) * t92 + qJ(3) * t37, t138 * (-pkin(7) * t72 + t186) + t140 * (-pkin(3) * t94 + pkin(7) * t73 + t191) - pkin(2) * t94 + qJ(3) * t45, t138 * (-pkin(7) * t65 - t18) + t140 * (-pkin(3) * t77 + pkin(7) * t66 + t19) - pkin(2) * t77 + qJ(3) * t40, -pkin(7) * t194 + t140 * (-pkin(3) * t71 + pkin(7) * t19) - pkin(2) * t71 + qJ(3) * t7, t138 * (t147 * t47 + t169) + t140 * (t144 * t47 - t168), t138 * (t144 * t76 + t147 * t30) + t140 * (t144 * t30 - t147 * t76), t138 * (t144 * t52 + t147 * t43) + t140 * (t144 * t43 - t147 * t52), t138 * (t147 * t46 - t169) + t140 * (t144 * t46 + t168), t138 * (-t144 * t48 + t147 * t44) + t140 * (t144 * t44 + t147 * t48), t138 * (t144 * t86 + t147 * t56) + t140 * (t144 * t56 - t147 * t86), t138 * (-pkin(7) * t22 - t144 * t9 + t147 * t16) + t140 * (-pkin(3) * t38 + pkin(7) * t23 + t144 * t16 + t147 * t9) - pkin(2) * t38 + qJ(3) * t11, t138 * (-pkin(7) * t24 - t10 * t144 + t147 * t17) + t140 * (-pkin(3) * t41 + pkin(7) * t25 + t10 * t147 + t144 * t17) - pkin(2) * t41 + qJ(3) * t12, t138 * (-pkin(7) * t20 + t147 * t4) + t140 * (pkin(7) * t21 + t144 * t4) + qJ(3) * t8 + (pkin(4) * t180 + t140 * t167 - pkin(2)) * t28, (t138 * (pkin(4) * t144 - pkin(8) * t147) + t140 * (-pkin(8) * t144 + t167) - pkin(2)) * t5 + t201 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, t171, -t120, t80, 0, 0, 0, 0, 0, 0, t92, t94, t77, t71, 0, 0, 0, 0, 0, 0, t38, t41, t28, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t111 - t110, t112, -t96, -t162, qJDD(4), -t32, -t33, 0, 0, t101 * t181 + t143 * t69, t143 * t50 + t146 * t51, t146 * t82 + t206, t146 * t68 + t99 * t182, t143 * t81 + t188, (-t101 * t146 - t143 * t99) * t109, pkin(4) * t50 + pkin(8) * t39 - t189, pkin(4) * t53 + pkin(8) * t42 + t193, pkin(4) * t67 + pkin(8) * t29 + t6, -pkin(4) * t26 + pkin(8) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t76, t52, -t78, -t48, t86, -t14, -t15, 0, 0;];
tauJ_reg = t13;
