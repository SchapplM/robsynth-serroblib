% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:44
% EndTime: 2019-12-05 16:19:47
% DurationCPUTime: 0.87s
% Computational Cost: add. (1252->188), mult. (2806->262), div. (0->0), fcn. (2079->10), ass. (0->116)
t102 = qJD(3) + qJD(5);
t108 = sin(qJ(5));
t110 = cos(qJ(5));
t105 = sin(pkin(9));
t106 = cos(pkin(9));
t109 = sin(qJ(3));
t111 = cos(qJ(3));
t74 = -t105 * t109 + t106 * t111;
t69 = t74 * qJD(2);
t61 = t110 * t69;
t75 = t105 * t111 + t106 * t109;
t71 = t75 * qJD(2);
t35 = -t108 * t71 + t61;
t142 = t35 * t102;
t139 = qJD(5) * t108;
t70 = t75 * qJD(3);
t40 = -qJD(2) * t70 + t74 * qJDD(2);
t138 = qJD(2) * qJD(3);
t134 = t111 * t138;
t135 = t109 * t138;
t41 = t75 * qJDD(2) - t105 * t135 + t106 * t134;
t5 = qJD(5) * t61 + t108 * t40 + t110 * t41 - t71 * t139;
t161 = t5 - t142;
t123 = t108 * t69 + t110 * t71;
t160 = t123 * t35;
t143 = t123 * t102;
t6 = t123 * qJD(5) + t108 * t41 - t110 * t40;
t159 = -t6 + t143;
t101 = pkin(8) + qJ(2);
t95 = sin(t101);
t96 = cos(t101);
t130 = g(1) * t96 + g(2) * t95;
t158 = t123 ^ 2 - t35 ^ 2;
t150 = t69 * pkin(7);
t146 = qJ(4) + pkin(6);
t87 = t146 * t111;
t66 = t109 * qJD(1) + qJD(2) * t87;
t144 = t106 * t66;
t145 = qJD(3) * pkin(3);
t86 = t146 * t109;
t64 = t111 * qJD(1) - qJD(2) * t86;
t60 = t64 + t145;
t24 = t105 * t60 + t144;
t13 = t24 + t150;
t94 = t111 * pkin(3) + pkin(2);
t83 = -t94 * qJD(2) + qJD(4);
t44 = -t69 * pkin(4) + t83;
t98 = qJ(3) + pkin(9) + qJ(5);
t91 = sin(t98);
t92 = cos(t98);
t157 = g(3) * t91 + t13 * t139 + t130 * t92 - t44 * t35;
t155 = qJD(5) - t102;
t132 = t146 * qJD(3);
t128 = qJD(2) * t132;
t153 = qJD(3) * qJD(1) + qJD(2) * qJD(4) + t146 * qJDD(2);
t97 = t111 * qJDD(1);
t22 = qJDD(3) * pkin(3) - t153 * t109 - t111 * t128 + t97;
t25 = (qJDD(1) - t128) * t109 + t153 * t111;
t7 = -t105 * t25 + t106 * t22;
t2 = qJDD(3) * pkin(4) - t41 * pkin(7) + t7;
t8 = t105 * t22 + t106 * t25;
t3 = t40 * pkin(7) + t8;
t154 = -g(3) * t92 - t108 * t3 + t110 * t2 - t44 * t123 + t130 * t91;
t149 = t71 * pkin(7);
t148 = pkin(3) * t105;
t147 = g(3) * t111;
t53 = t105 * t66;
t28 = t106 * t64 - t53;
t65 = t111 * qJD(4) - t109 * t132;
t67 = -t109 * qJD(4) - t111 * t132;
t29 = t105 * t67 + t106 * t65;
t46 = -t105 * t86 + t106 * t87;
t141 = qJDD(1) - g(3);
t103 = t109 ^ 2;
t140 = -t111 ^ 2 + t103;
t137 = t111 * qJDD(2);
t136 = t109 * t145;
t23 = t106 * t60 - t53;
t26 = -t105 * t64 - t144;
t27 = -t105 * t65 + t106 * t67;
t45 = -t105 * t87 - t106 * t86;
t129 = g(1) * t95 - g(2) * t96;
t100 = qJDD(3) + qJDD(5);
t43 = t108 * t74 + t110 * t75;
t42 = t108 * t75 - t110 * t74;
t73 = t74 * qJD(3);
t9 = -t42 * qJD(5) - t108 * t70 + t110 * t73;
t127 = t43 * t100 + t9 * t102;
t12 = qJD(3) * pkin(4) - t149 + t23;
t126 = -t108 * t12 - t110 * t13;
t30 = -t75 * pkin(7) + t45;
t31 = t74 * pkin(7) + t46;
t125 = -t108 * t31 + t110 * t30;
t124 = t108 * t30 + t110 * t31;
t93 = t106 * pkin(3) + pkin(4);
t121 = t108 * t93 + t110 * t148;
t120 = -t108 * t148 + t110 * t93;
t119 = -0.2e1 * pkin(2) * t138 - pkin(6) * qJDD(3);
t117 = pkin(3) * t135 - t94 * qJDD(2) + qJDD(4);
t112 = qJD(3) ^ 2;
t116 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t112 + t129;
t113 = qJD(2) ^ 2;
t115 = t113 * pkin(2) - qJDD(2) * pkin(6) + t130;
t85 = qJDD(3) * t111 - t112 * t109;
t84 = qJDD(3) * t109 + t112 * t111;
t52 = -t74 * pkin(4) - t94;
t48 = t70 * pkin(4) + t136;
t47 = t109 * qJD(2) * pkin(3) + t71 * pkin(4);
t18 = -t40 * pkin(4) + t117;
t17 = -t70 * pkin(7) + t29;
t16 = t28 - t149;
t15 = -t73 * pkin(7) + t27;
t14 = t26 - t150;
t10 = t43 * qJD(5) + t108 * t73 + t110 * t70;
t4 = -t10 * t102 - t42 * t100;
t1 = [t141, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t84, t75 * t40 - t74 * t41 + t73 * t69 + t70 * t71, -t23 * t70 + t24 * t73 + t7 * t74 + t8 * t75 - g(3), 0, 0, 0, 0, 0, t4, -t127; 0, qJDD(2), t129, t130, t103 * qJDD(2) + 0.2e1 * t109 * t134, 0.2e1 * t109 * t137 - 0.2e1 * t140 * t138, t84, t85, 0, t119 * t109 + t116 * t111, -t116 * t109 + t119 * t111, -t23 * t73 - t24 * t70 - t27 * t71 + t29 * t69 + t46 * t40 - t45 * t41 - t7 * t75 + t8 * t74 - t130, t8 * t46 + t24 * t29 + t7 * t45 + t23 * t27 - t117 * t94 + t83 * t136 - g(1) * (t146 * t96 - t95 * t94) - g(2) * (t146 * t95 + t96 * t94), t123 * t9 + t5 * t43, -t10 * t123 + t35 * t9 - t5 * t42 - t43 * t6, t127, t4, 0, -t48 * t35 + t52 * t6 + t18 * t42 + t44 * t10 + (-t124 * qJD(5) - t108 * t17 + t110 * t15) * t102 + t125 * t100 + t129 * t92, t48 * t123 + t52 * t5 + t18 * t43 + t44 * t9 - (t125 * qJD(5) + t108 * t15 + t110 * t17) * t102 - t124 * t100 - t129 * t91; 0, 0, 0, 0, -t109 * t113 * t111, t140 * t113, t109 * qJDD(2), t137, qJDD(3), t109 * t115 - t147 + t97, -t141 * t109 + t115 * t111, (t24 + t26) * t71 + (t23 - t28) * t69 + (t105 * t40 - t106 * t41) * pkin(3), -t23 * t26 - t24 * t28 + (-t147 + t105 * t8 + t106 * t7 + (-qJD(2) * t83 + t130) * t109) * pkin(3), -t160, t158, t161, t159, t100, t120 * t100 + t47 * t35 - (-t108 * t16 + t110 * t14) * t102 + (-t121 * t102 + t126) * qJD(5) + t154, -t121 * t100 - t110 * t3 - t108 * t2 - t47 * t123 + (t108 * t14 + t110 * t16) * t102 + (-t120 * t102 - t110 * t12) * qJD(5) + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 ^ 2 - t71 ^ 2, t23 * t71 - t24 * t69 + t117 - t129, 0, 0, 0, 0, 0, t6 + t143, t5 + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, t158, t161, t159, t100, t155 * t126 + t154, (-t13 * t102 - t2) * t108 + (-t155 * t12 - t3) * t110 + t157;];
tau_reg = t1;
