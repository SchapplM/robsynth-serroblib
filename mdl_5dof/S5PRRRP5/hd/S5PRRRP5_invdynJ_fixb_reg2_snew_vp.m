% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRP5
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:18
% EndTime: 2019-12-05 16:49:24
% DurationCPUTime: 1.43s
% Computational Cost: add. (3350->207), mult. (7018->269), div. (0->0), fcn. (4727->8), ass. (0->136)
t129 = qJDD(3) + qJDD(4);
t137 = sin(qJ(4));
t140 = cos(qJ(4));
t141 = cos(qJ(3));
t138 = sin(qJ(3));
t160 = qJD(2) * t138;
t102 = -t140 * t141 * qJD(2) + t137 * t160;
t104 = (t137 * t141 + t138 * t140) * qJD(2);
t79 = t104 * t102;
t184 = t79 - t129;
t185 = t184 * pkin(4);
t183 = t137 * t184;
t182 = t140 * t184;
t157 = qJD(2) * qJD(3);
t151 = t141 * t157;
t156 = t138 * qJDD(2);
t109 = t151 + t156;
t125 = t141 * qJDD(2);
t152 = t138 * t157;
t110 = t125 - t152;
t67 = -qJD(4) * t102 + t109 * t140 + t110 * t137;
t130 = qJD(3) + qJD(4);
t94 = t130 * t102;
t181 = t67 - t94;
t144 = qJD(2) ^ 2;
t119 = t138 * t144 * t141;
t116 = qJDD(3) + t119;
t134 = sin(pkin(8));
t135 = cos(pkin(8));
t114 = -g(1) * t134 + g(2) * t135;
t115 = -g(1) * t135 - g(2) * t134;
t139 = sin(qJ(2));
t142 = cos(qJ(2));
t161 = -g(3) + qJDD(1);
t92 = t142 * t115 + t139 * t161;
t82 = -t144 * pkin(2) + qJDD(2) * pkin(6) + t92;
t70 = -t141 * t114 + t138 * t82;
t47 = (-t109 + t151) * pkin(7) + t116 * pkin(3) - t70;
t118 = qJD(3) * pkin(3) - pkin(7) * t160;
t132 = t141 ^ 2;
t127 = t132 * t144;
t71 = t114 * t138 + t141 * t82;
t48 = -pkin(3) * t127 + pkin(7) * t110 - qJD(3) * t118 + t71;
t25 = t137 * t48 - t140 * t47;
t179 = qJ(5) * t94 + 0.2e1 * qJD(5) * t104 + t185 + t25;
t100 = t102 ^ 2;
t28 = t137 * t47 + t140 * t48;
t150 = t137 * t109 - t140 * t110;
t66 = -qJD(4) * t104 - t150;
t87 = pkin(4) * t130 - qJ(5) * t104;
t17 = -t100 * pkin(4) + t66 * qJ(5) - 0.2e1 * qJD(5) * t102 - t130 * t87 + t28;
t101 = t104 ^ 2;
t128 = t130 ^ 2;
t146 = (-qJD(4) + t130) * t104 - t150;
t59 = t67 + t94;
t32 = t137 * t146 - t140 * t59;
t178 = pkin(7) * t32;
t73 = -t128 - t100;
t42 = t137 * t73 - t182;
t177 = pkin(7) * t42;
t75 = t79 + t129;
t170 = t137 * t75;
t86 = -t101 - t128;
t60 = t140 * t86 - t170;
t176 = pkin(7) * t60;
t9 = t137 * t28 - t140 * t25;
t175 = t138 * t9;
t33 = t137 * t59 + t140 * t146;
t13 = -t138 * t32 + t141 * t33;
t68 = -t100 - t101;
t174 = -pkin(2) * t68 + pkin(6) * t13;
t43 = t140 * t73 + t183;
t23 = -t138 * t42 + t141 * t43;
t54 = (qJD(4) + t130) * t104 + t150;
t173 = -pkin(2) * t54 + pkin(6) * t23;
t168 = t140 * t75;
t61 = -t137 * t86 - t168;
t34 = -t138 * t60 + t141 * t61;
t172 = -pkin(2) * t181 + pkin(6) * t34;
t91 = -t115 * t139 + t142 * t161;
t81 = -qJDD(2) * pkin(2) - t144 * pkin(6) - t91;
t62 = -t110 * pkin(3) - pkin(7) * t127 + t118 * t160 + t81;
t171 = t137 * t62;
t169 = t140 * t62;
t167 = qJ(5) * t137;
t166 = qJ(5) * t140;
t165 = t130 * t137;
t164 = t130 * t140;
t163 = t138 * t116;
t162 = t141 * (qJDD(3) - t119);
t155 = -pkin(3) * t68 + pkin(7) * t33;
t154 = -pkin(3) * t54 + pkin(7) * t43;
t153 = -pkin(3) * t181 + pkin(7) * t61;
t10 = t137 * t25 + t140 * t28;
t39 = t138 * t70 + t141 * t71;
t148 = pkin(4) * t86 - t17;
t15 = -qJ(5) * t67 - t179;
t145 = t15 - t185;
t21 = -t66 * pkin(4) - t100 * qJ(5) + t104 * t87 + qJDD(5) + t62;
t143 = qJD(3) ^ 2;
t131 = t138 ^ 2;
t126 = t131 * t144;
t113 = t126 + t127;
t112 = (t131 + t132) * qJDD(2);
t111 = t125 - 0.2e1 * t152;
t108 = 0.2e1 * t151 + t156;
t89 = -t101 + t128;
t88 = t100 - t128;
t85 = -t162 - t138 * (-t126 - t143);
t84 = t141 * (-t127 - t143) - t163;
t77 = t101 - t100;
t53 = pkin(3) * t60;
t49 = pkin(4) * t59;
t41 = pkin(3) * t42;
t38 = (t138 * (-t102 * t140 + t104 * t137) + t141 * (-t102 * t137 - t104 * t140)) * t130;
t37 = -pkin(4) * t181 - qJ(5) * t75;
t36 = t138 * (t140 * t88 - t170) + t141 * (t137 * t88 + t168);
t35 = t138 * (-t137 * t89 - t182) + t141 * (t140 * t89 - t183);
t30 = pkin(3) * t32;
t27 = t138 * (-t104 * t165 + t140 * t67) + t141 * (t104 * t164 + t137 * t67);
t26 = t138 * (t102 * t164 - t137 * t66) + t141 * (t102 * t165 + t140 * t66);
t20 = -qJ(5) * t86 + t21;
t19 = t139 * t34 - t142 * t181;
t18 = t139 * t23 - t142 * t54;
t16 = -pkin(4) * t54 + qJ(5) * t73 - t21;
t14 = pkin(4) * t15;
t12 = t138 * (-t137 * t181 - t140 * t54) + t141 * (-t137 * t54 + t140 * t181);
t8 = t13 * t139 - t142 * t68;
t7 = (t59 + t67) * qJ(5) + t179;
t6 = -pkin(4) * t68 + qJ(5) * t146 + t17;
t5 = -pkin(4) * t21 + qJ(5) * t17;
t4 = -t137 * t15 + t140 * t17;
t3 = t137 * t17 + t140 * t15;
t2 = t10 * t141 - t175;
t1 = -t138 * t3 + t141 * t4;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t161, 0, 0, 0, 0, 0, 0, qJDD(2) * t142 - t139 * t144, -qJDD(2) * t139 - t142 * t144, 0, t139 * t92 + t142 * t91, 0, 0, 0, 0, 0, 0, t111 * t142 + t139 * t84, -t108 * t142 + t139 * t85, t112 * t139 + t113 * t142, t139 * t39 - t142 * t81, 0, 0, 0, 0, 0, 0, t18, t19, t8, t139 * t2 - t142 * t62, 0, 0, 0, 0, 0, 0, t18, t19, t8, t1 * t139 - t142 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t91, -t92, 0, 0, (t109 + t151) * t138, t108 * t141 + t111 * t138, t163 + t141 * (-t126 + t143), (t110 - t152) * t141, t138 * (t127 - t143) + t162, 0, pkin(2) * t111 + pkin(6) * t84 - t141 * t81, -pkin(2) * t108 + pkin(6) * t85 + t138 * t81, pkin(2) * t113 + pkin(6) * t112 + t39, -pkin(2) * t81 + pkin(6) * t39, t27, t12, t35, t26, t36, t38, t138 * (t171 - t177) + t141 * (t154 - t169) + t173, t138 * (t169 - t176) + t141 * (t153 + t171) + t172, t138 * (-t9 - t178) + t141 * (t10 + t155) + t174, -pkin(7) * t175 + t141 * (-pkin(3) * t62 + pkin(7) * t10) - pkin(2) * t62 + pkin(6) * t2, t27, t12, t35, t26, t36, t38, t138 * (-t137 * t16 + t166 * t184 - t177) + t141 * (t140 * t16 + t167 * t184 + t154) + t173, t138 * (-t137 * t37 + t140 * t20 - t176) + t141 * (t137 * t20 + t140 * t37 + t153) + t172, t138 * (-t137 * t6 + t140 * t7 - t178) + t141 * (t137 * t7 + t140 * t6 + t155) + t174, t138 * (-pkin(7) * t3 - t137 * t5 - t15 * t166) + t141 * (-pkin(3) * t21 + pkin(7) * t4 + t140 * t5 - t15 * t167) - pkin(2) * t21 + pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t126 - t127, t156, t119, t125, qJDD(3), -t70, -t71, 0, 0, t79, t77, t59, -t79, t146, t129, -t25 + t41, t53 - t28, t30, pkin(3) * t9, t79, t77, t59, -t79, t146, t129, t145 + t41, t53 + t148, -t49 + t30, pkin(3) * t3 + t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t77, t59, -t79, t146, t129, -t25, -t28, 0, 0, t79, t77, t59, -t79, t146, t129, t145, t148, -t49, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t181, t68, t21;];
tauJ_reg = t11;
