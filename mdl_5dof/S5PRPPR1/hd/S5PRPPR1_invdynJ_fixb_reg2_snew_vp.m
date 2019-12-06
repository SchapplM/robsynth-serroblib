% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPPR1
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:10
% EndTime: 2019-12-05 15:22:17
% DurationCPUTime: 1.74s
% Computational Cost: add. (3326->212), mult. (8390->327), div. (0->0), fcn. (5776->10), ass. (0->150)
t109 = sin(pkin(8));
t103 = t109 ^ 2;
t112 = cos(pkin(8));
t105 = t112 ^ 2;
t182 = t103 + t105;
t131 = -t112 * pkin(3) - t109 * qJ(4);
t152 = qJD(2) * t109;
t118 = qJD(2) ^ 2;
t106 = t118 * qJ(3);
t107 = qJDD(2) * pkin(2);
t115 = sin(qJ(2));
t117 = cos(qJ(2));
t110 = sin(pkin(7));
t113 = cos(pkin(7));
t88 = g(1) * t110 - g(2) * t113;
t89 = -g(1) * t113 - g(2) * t110;
t136 = -t115 * t89 + t117 * t88;
t59 = qJDD(3) - t106 - t107 - t136;
t181 = -0.2e1 * qJD(4) * t152 + t131 * qJDD(2) + t59;
t108 = sin(pkin(9));
t111 = cos(pkin(9));
t114 = sin(qJ(5));
t116 = cos(qJ(5));
t128 = t108 * t116 + t111 * t114;
t123 = t128 * t109;
t67 = qJD(2) * t123;
t151 = qJD(2) * t112;
t94 = -qJD(5) + t151;
t172 = t67 * t94;
t148 = t109 * qJDD(2);
t138 = t108 * t148;
t157 = t109 * t111;
t140 = t116 * t157;
t51 = -qJD(5) * t67 + qJDD(2) * t140 - t114 * t138;
t180 = t51 + t172;
t174 = 2 * qJD(3);
t130 = -t115 * t88 - t117 * t89;
t60 = -pkin(2) * t118 + qJDD(2) * qJ(3) - t130;
t179 = qJD(2) * t174 + t60;
t69 = -t114 * t108 * t152 + qJD(2) * t140;
t178 = (qJD(5) + t94) * t69;
t171 = t69 * t67;
t147 = t112 * qJDD(2);
t93 = -qJDD(5) + t147;
t124 = -t93 - t171;
t177 = t114 * t124;
t176 = t116 * t124;
t175 = t182 * t106 - t107 + t59;
t65 = t67 ^ 2;
t66 = t69 ^ 2;
t92 = t94 ^ 2;
t125 = t181 * t111;
t127 = -t112 * pkin(4) - pkin(6) * t157;
t161 = t103 * t111;
t155 = -g(3) + qJDD(1);
t46 = t109 * t155 + t179 * t112;
t81 = t131 * qJD(2);
t36 = t81 * t151 + t46;
t17 = t127 * qJDD(2) + (-t36 + (pkin(6) * t109 * t112 - pkin(4) * t161) * t118) * t108 + t125;
t22 = t181 * t108 + t111 * t36;
t77 = t127 * qJD(2);
t102 = t108 ^ 2;
t160 = t103 * t118;
t91 = t102 * t160;
t18 = -pkin(4) * t91 - pkin(6) * t138 + t77 * t151 + t22;
t6 = t114 * t18 - t116 * t17;
t7 = t114 * t17 + t116 * t18;
t3 = t114 * t7 - t116 * t6;
t173 = t111 * t3;
t158 = t108 * t118;
t133 = t158 * t161;
t74 = -t133 + t147;
t170 = t108 * t74;
t75 = -t133 - t147;
t169 = t111 * t75;
t99 = t112 * t155;
t150 = qJDD(4) - t99;
t154 = t174 + t81;
t25 = -pkin(6) * t91 + (pkin(4) * qJDD(2) * t108 + t60 + (t111 * t77 + t154) * qJD(2)) * t109 + t150;
t168 = t114 * t25;
t40 = t93 - t171;
t167 = t114 * t40;
t166 = t114 * t94;
t165 = t116 * t25;
t164 = t116 * t40;
t163 = t116 * t94;
t104 = t111 ^ 2;
t162 = t103 * t104;
t159 = t105 * t118;
t156 = t112 * t118;
t149 = qJDD(2) * t111;
t146 = t112 * t171;
t142 = t104 * t160;
t141 = t108 * t156;
t139 = t111 * t156;
t4 = t114 * t6 + t116 * t7;
t45 = t179 * t109 - t99;
t137 = t109 * t45 + t112 * t46;
t132 = t108 * t139;
t21 = t108 * t36 - t125;
t10 = t108 * t22 - t111 * t21;
t129 = t103 * t132;
t126 = -pkin(2) + t131;
t122 = t128 * t148;
t121 = qJDD(2) * t123;
t101 = t105 * qJDD(2);
t100 = t103 * qJDD(2);
t83 = t182 * t118;
t82 = t109 * t139;
t79 = (-t105 - t162) * t118;
t78 = -t91 - t159;
t76 = t91 + t142;
t73 = (t141 - t149) * t109;
t72 = (t141 + t149) * t109;
t71 = -t82 + t138;
t70 = t82 + t138;
t62 = -t66 + t92;
t61 = t65 - t92;
t57 = -t66 - t92;
t56 = -t108 * t79 + t111 * t74;
t55 = -t108 * t75 + t111 * t78;
t54 = t111 * t79 + t170;
t53 = t108 * t78 + t169;
t52 = t66 - t65;
t50 = -qJD(5) * t69 - t121;
t48 = -t108 * t73 - t111 * t70;
t47 = -t108 * t70 + t111 * t73;
t39 = -t92 - t65;
t37 = -t65 - t66;
t35 = (t154 * qJD(2) + t60) * t109 + t150;
t33 = t51 - t172;
t30 = -t121 - t178;
t29 = t122 + t178;
t28 = (qJD(5) - t94) * t69 + t122;
t27 = -t114 * t57 + t164;
t26 = t116 * t57 + t167;
t24 = t116 * t39 - t177;
t23 = t114 * t39 + t176;
t20 = t114 * t33 + t116 * t30;
t19 = t114 * t30 - t116 * t33;
t15 = -t108 * t26 + t111 * t27;
t14 = t108 * t27 + t111 * t26;
t13 = -t108 * t23 + t111 * t24;
t12 = t108 * t24 + t111 * t23;
t11 = t108 * t21 + t111 * t22;
t9 = -t108 * t19 + t111 * t20;
t8 = t108 * t20 + t111 * t19;
t2 = -t108 * t3 + t111 * t4;
t1 = t108 * t4 + t173;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t155, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t46 - t112 * t45, 0, 0, 0, 0, 0, 0, t109 * t55 - t112 * t71, t109 * t56 - t112 * t72, t109 * t48 + t112 * t76, t109 * t11 - t112 * t35, 0, 0, 0, 0, 0, 0, t109 * t13 - t112 * t28, t109 * t15 - t112 * t180, t109 * t9 - t112 * t37, t109 * t2 - t112 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t136, t130, 0, 0, t100, 0.2e1 * t109 * t147, 0, t101, 0, 0, -t175 * t112, t175 * t109, pkin(2) * t83 + qJ(3) * (t101 + t100) + t137, -pkin(2) * t59 + qJ(3) * t137, -t129 + (qJDD(2) * t104 + t132) * t103, t109 * (-t108 * t72 - t111 * t71) + t112 * (t91 - t142), t109 * (t169 - (t105 - t162) * t158) + t112 * t73, t129 + (qJDD(2) * t102 - t132) * t103, t109 * (t111 * (t91 - t159) + t170) + t112 * t70, t101, t109 * (-qJ(4) * t53 + t108 * t35) + t112 * (-pkin(3) * t53 + t21) - pkin(2) * t53 + qJ(3) * (t109 * t71 + t112 * t55), t109 * (-qJ(4) * t54 + t111 * t35) + t112 * (-pkin(3) * t54 + t22) - pkin(2) * t54 + qJ(3) * (t109 * t72 + t112 * t56), -t109 * t10 + qJ(3) * (-t109 * t76 + t112 * t48) + t126 * t47, qJ(3) * (t109 * t35 + t11 * t112) + t126 * t10, t109 * (t111 * (t116 * t51 + t69 * t166) - t108 * (t114 * t51 - t69 * t163)) - t146, t109 * (t111 * (-t114 * t180 - t116 * t28) - t108 * (-t114 * t28 + t116 * t180)) - t112 * t52, t109 * (t111 * (-t114 * t62 + t176) - t108 * (t116 * t62 + t177)) - t112 * t33, t109 * (t111 * (-t114 * t50 - t67 * t163) - t108 * (t116 * t50 - t67 * t166)) + t146, t109 * (t111 * (t116 * t61 + t167) - t108 * (t114 * t61 - t164)) + t112 * t29, t112 * t93 + t109 * (t111 * (-t114 * t69 + t116 * t67) - t108 * (t114 * t67 + t116 * t69)) * t94, t109 * (t111 * (-pkin(6) * t23 + t168) - t108 * (-pkin(4) * t28 + pkin(6) * t24 - t165) - qJ(4) * t12) + t112 * (-pkin(3) * t12 - pkin(4) * t23 + t6) - pkin(2) * t12 + qJ(3) * (t109 * t28 + t112 * t13), t109 * (t111 * (-pkin(6) * t26 + t165) - t108 * (-pkin(4) * t180 + pkin(6) * t27 + t168) - qJ(4) * t14) + t112 * (-pkin(3) * t14 - pkin(4) * t26 + t7) - pkin(2) * t14 + qJ(3) * (t109 * t180 + t112 * t15), t109 * (t111 * (-pkin(6) * t19 - t3) - t108 * (-pkin(4) * t37 + pkin(6) * t20 + t4) - qJ(4) * t8) + t112 * (-pkin(3) * t8 - pkin(4) * t19) - pkin(2) * t8 + qJ(3) * (t109 * t37 + t112 * t9), t109 * (-pkin(6) * t173 - t108 * (-pkin(4) * t25 + pkin(6) * t4) - qJ(4) * t1) + t112 * (-pkin(3) * t1 - pkin(4) * t3) - pkin(2) * t1 + qJ(3) * (t109 * t25 + t112 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t148, -t83, t59, 0, 0, 0, 0, 0, 0, t53, t54, t47, t10, 0, 0, 0, 0, 0, 0, t12, t14, t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t72, -t76, t35, 0, 0, 0, 0, 0, 0, t28, t180, t37, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t52, t33, -t171, -t29, -t93, -t6, -t7, 0, 0;];
tauJ_reg = t5;
