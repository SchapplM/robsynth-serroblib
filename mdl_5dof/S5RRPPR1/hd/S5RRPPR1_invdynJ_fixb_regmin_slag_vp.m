% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:30
% EndTime: 2019-12-05 18:18:34
% DurationCPUTime: 0.89s
% Computational Cost: add. (1086->178), mult. (1749->230), div. (0->0), fcn. (1162->16), ass. (0->124)
t112 = cos(qJ(2));
t155 = pkin(1) * qJD(2);
t141 = qJD(1) * t155;
t109 = sin(qJ(2));
t146 = qJDD(1) * t109;
t172 = pkin(1) * t146 + t112 * t141;
t106 = cos(pkin(9));
t111 = cos(qJ(5));
t148 = t111 * t106;
t104 = sin(pkin(9));
t108 = sin(qJ(5));
t151 = t108 * t104;
t52 = -t148 + t151;
t53 = t111 * t104 + t108 * t106;
t103 = qJ(1) + qJ(2);
t93 = pkin(8) + t103;
t79 = sin(t93);
t80 = cos(t93);
t171 = g(2) * t80 + g(3) * t79;
t94 = sin(t103);
t95 = cos(t103);
t170 = g(2) * t95 + g(3) * t94;
t102 = qJD(1) + qJD(2);
t39 = t53 * t102;
t153 = t104 ^ 2 + t106 ^ 2;
t129 = -g(2) * t79 + g(3) * t80;
t107 = cos(pkin(8));
t156 = pkin(1) * qJD(1);
t142 = t112 * t156;
t105 = sin(pkin(8));
t143 = t109 * t156;
t70 = t105 * t143;
t46 = t107 * t142 - t70;
t147 = qJD(4) - t46;
t169 = pkin(2) * t94;
t168 = pkin(2) * t95;
t164 = pkin(1) * t109;
t163 = t106 * pkin(4);
t162 = t107 * pkin(2);
t110 = sin(qJ(1));
t161 = t110 * pkin(1);
t160 = t112 * pkin(1);
t113 = cos(qJ(1));
t159 = t113 * pkin(1);
t86 = qJDD(1) * t160;
t98 = qJDD(1) + qJDD(2);
t41 = t98 * pkin(2) - t109 * t141 + t86;
t22 = t105 * t41 + t172 * t107;
t14 = t98 * qJ(4) + t102 * qJD(4) + t22;
t9 = t104 * qJDD(3) + t106 * t14;
t158 = t171 * t106;
t56 = t102 * pkin(2) + t142;
t71 = t107 * t143;
t34 = t105 * t56 + t71;
t152 = t107 * t109;
t85 = pkin(2) + t160;
t157 = pkin(1) * t152 + t105 * t85;
t154 = t106 * t98;
t138 = t102 * t148;
t145 = qJD(5) * t138 + t53 * t98;
t144 = t86 + t170;
t75 = t105 * t164;
t140 = t102 * t151;
t137 = -pkin(3) - t163;
t136 = -g(2) * t94 + g(3) * t95;
t21 = -t172 * t105 + t107 * t41;
t117 = qJDD(4) - t21;
t10 = t137 * t98 + t117;
t33 = t107 * t56 - t70;
t127 = qJD(4) - t33;
t25 = t102 * t137 + t127;
t48 = t53 * qJD(5);
t101 = pkin(9) + qJ(5);
t92 = cos(t101);
t135 = t10 * t52 + t171 * t92 + t25 * t48;
t134 = t107 * t85 - t75;
t133 = qJD(1) * (-qJD(2) + t102);
t132 = qJD(2) * (-qJD(1) - t102);
t43 = -pkin(3) - t134;
t128 = t52 * t98;
t45 = (t105 * t112 + t152) * t155;
t125 = -t102 * t45 - t43 * t98;
t44 = t105 * t142 + t71;
t81 = -pkin(3) - t162;
t124 = t102 * t44 - t81 * t98;
t42 = qJ(4) + t157;
t31 = (-pkin(7) - t42) * t104;
t96 = t106 * pkin(7);
t32 = t106 * t42 + t96;
t123 = t108 * t32 - t111 * t31;
t122 = t108 * t31 + t111 * t32;
t77 = t105 * pkin(2) + qJ(4);
t49 = (-pkin(7) - t77) * t104;
t50 = t106 * t77 + t96;
t121 = t108 * t50 - t111 * t49;
t120 = t108 * t49 + t111 * t50;
t119 = -t79 * pkin(3) + t80 * qJ(4) - t169;
t118 = t107 * t112 * t155 - qJD(2) * t75;
t88 = t106 * qJDD(3);
t8 = -t104 * t14 + t88;
t116 = -t8 * t104 + t9 * t106 - t129;
t115 = -t80 * pkin(3) - t79 * qJ(4) - t168;
t16 = -t98 * pkin(3) + t117;
t47 = t52 * qJD(5);
t91 = sin(t101);
t114 = t10 * t53 - t171 * t91 - t25 * t47;
t58 = t137 - t162;
t40 = qJD(4) + t118;
t37 = -t138 + t140;
t35 = t43 - t163;
t30 = t102 * qJ(4) + t34;
t29 = -t102 * pkin(3) + t127;
t27 = -t48 * qJD(5) - t52 * qJDD(5);
t26 = -t47 * qJD(5) + t53 * qJDD(5);
t24 = t104 * qJD(3) + t106 * t30;
t23 = t106 * qJD(3) - t104 * t30;
t20 = t102 * t48 + t128;
t19 = -qJD(5) * t140 + t145;
t15 = t16 * t104;
t4 = pkin(7) * t154 + t9;
t3 = t88 + (-pkin(7) * t98 - t14) * t104;
t2 = t19 * t53 - t39 * t47;
t1 = -t19 * t52 - t53 * t20 + t47 * t37 - t39 * t48;
t5 = [qJDD(1), g(2) * t113 + g(3) * t110, -g(2) * t110 + g(3) * t113, t98, (t109 * t132 + t112 * t98) * pkin(1) + t144, ((-qJDD(1) - t98) * t109 + t112 * t132) * pkin(1) + t136, t22 * t157 + t34 * t118 + t21 * t134 - t33 * t45 - g(2) * (-t159 - t168) - g(3) * (-t161 - t169), (t125 - t16) * t106 + t158, t15 + (-t125 - t171) * t104, t116 + t153 * (t102 * t40 + t42 * t98), t16 * t43 + t29 * t45 - g(2) * (t115 - t159) - g(3) * (t119 - t161) + (t24 * t40 + t9 * t42) * t106 + (-t23 * t40 - t8 * t42) * t104, t2, t1, t26, t27, 0, t45 * t37 + t35 * t20 - t123 * qJDD(5) + (-qJD(5) * t122 - t40 * t53) * qJD(5) + t135, t45 * t39 + t35 * t19 - t122 * qJDD(5) + (qJD(5) * t123 + t40 * t52) * qJD(5) + t114; 0, 0, 0, t98, t133 * t164 + t144, (t112 * t133 - t146) * pkin(1) + t136, t33 * t44 - t34 * t46 + (t105 * t22 + t107 * t21 + t170) * pkin(2), (t124 - t16) * t106 + t158, t15 + (-t124 - t171) * t104, t116 + (t147 * t102 + t98 * t77) * t153, t16 * t81 - t29 * t44 - g(2) * t115 - g(3) * t119 + (t147 * t24 + t9 * t77) * t106 + (-t147 * t23 - t8 * t77) * t104, t2, t1, t26, t27, 0, t58 * t20 - t121 * qJDD(5) - t44 * t37 + (-qJD(5) * t120 - t147 * t53) * qJD(5) + t135, t58 * t19 - t120 * qJDD(5) - t44 * t39 + (qJD(5) * t121 + t147 * t52) * qJD(5) + t114; 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, t9 * t104 + t8 * t106 - g(1), 0, 0, 0, 0, 0, t27, -t26; 0, 0, 0, 0, 0, 0, 0, -t154, t104 * t98, -t153 * t102 ^ 2, (t23 * t104 - t24 * t106) * t102 + t16 - t171, 0, 0, 0, 0, 0, 0.2e1 * t39 * qJD(5) + t128, (-t37 - t140) * qJD(5) + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t37, -t37 ^ 2 + t39 ^ 2, (t37 - t140) * qJD(5) + t145, -t128, qJDD(5), -g(1) * t92 - t108 * t4 + t111 * t3 + t129 * t91 - t25 * t39, g(1) * t91 - t108 * t3 - t111 * t4 + t129 * t92 + t25 * t37;];
tau_reg = t5;
