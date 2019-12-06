% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:08:01
% EndTime: 2019-12-05 17:08:03
% DurationCPUTime: 0.90s
% Computational Cost: add. (1028->171), mult. (1507->237), div. (0->0), fcn. (1008->12), ass. (0->130)
t172 = pkin(7) + pkin(8);
t100 = cos(qJ(4));
t96 = sin(qJ(5));
t97 = sin(qJ(4));
t99 = cos(qJ(5));
t40 = t96 * t100 + t99 * t97;
t92 = qJD(2) + qJD(3);
t30 = t40 * t92;
t90 = pkin(9) + qJ(2);
t82 = qJ(3) + t90;
t72 = cos(t82);
t169 = g(2) * t72;
t150 = pkin(2) * qJD(2);
t98 = sin(qJ(3));
t132 = t98 * t150;
t101 = cos(qJ(3));
t165 = t101 * pkin(2);
t152 = -qJD(3) * t132 + qJDD(2) * t165;
t89 = qJDD(2) + qJDD(3);
t168 = t89 * pkin(3);
t32 = -t152 - t168;
t171 = t32 + t169;
t71 = sin(t82);
t153 = g(1) * t72 + g(2) * t71;
t125 = t172 * t92 + t132;
t23 = t100 * qJD(1) - t125 * t97;
t91 = qJD(4) + qJD(5);
t24 = t97 * qJD(1) + t100 * t125;
t68 = g(1) * t71;
t167 = t92 * pkin(3);
t75 = pkin(2) * t98 + pkin(7);
t166 = -pkin(8) - t75;
t145 = t99 * t100;
t131 = t92 * t145;
t158 = t96 * t97;
t135 = t92 * t158;
t28 = -t131 + t135;
t164 = t30 * t28;
t95 = qJ(4) + qJ(5);
t84 = sin(t95);
t163 = t71 * t84;
t85 = cos(t95);
t162 = t71 * t85;
t161 = t72 * t84;
t160 = t72 * t85;
t159 = t92 * t97;
t157 = t97 * t89;
t156 = t99 * t24;
t141 = t97 * qJD(4);
t139 = qJD(2) * t101;
t128 = pkin(2) * t139;
t51 = -t128 - t167;
t154 = t100 * t68 + t141 * t51;
t93 = t97 ^ 2;
t151 = -t100 ^ 2 + t93;
t149 = t100 * t89;
t148 = t100 * t97;
t147 = t101 * t91;
t144 = qJD(3) * t98;
t143 = qJD(5) * t96;
t140 = qJDD(1) - g(3);
t138 = qJD(3) * t101;
t137 = qJDD(2) * t98;
t136 = t100 * qJD(4);
t134 = t51 * t136 + t171 * t97;
t133 = pkin(4) * t141;
t130 = pkin(2) * t138;
t129 = t92 * t144;
t127 = t92 * t136;
t77 = -pkin(4) * t100 - pkin(3);
t33 = t89 * pkin(7) + (qJD(2) * t138 + t137) * pkin(2);
t126 = pkin(8) * t89 + t33;
t123 = qJD(4) * t172;
t122 = qJD(4) * t166;
t121 = t92 * t132;
t120 = -t145 * t89 + t157 * t96;
t19 = qJD(4) * pkin(4) + t23;
t118 = -t19 * t96 - t156;
t113 = t145 - t158;
t20 = t91 * t113;
t88 = qJDD(4) + qJDD(5);
t11 = t20 * t91 + t40 * t88;
t37 = t166 * t97;
t86 = t100 * pkin(8);
t38 = t100 * t75 + t86;
t117 = t37 * t99 - t38 * t96;
t116 = t37 * t96 + t38 * t99;
t65 = t172 * t97;
t66 = pkin(7) * t100 + t86;
t115 = -t65 * t99 - t66 * t96;
t114 = -t65 * t96 + t66 * t99;
t112 = t152 + t68 - t169;
t17 = (t141 * t92 - t149) * pkin(4) + t32;
t31 = t77 * t92 - t128;
t111 = -g(1) * t163 + g(2) * t161 + t17 * t40 + t20 * t31;
t21 = t91 * t40;
t110 = g(1) * t162 - g(2) * t160 - t113 * t17 + t21 * t31;
t109 = -t51 * t92 + t153 - t33;
t102 = qJD(4) ^ 2;
t108 = pkin(7) * t102 - t121 - t168;
t76 = -pkin(3) - t165;
t107 = pkin(2) * t129 + t102 * t75 + t76 * t89;
t8 = qJD(5) * t131 + t99 * t127 - t91 * t135 + t40 * t89;
t106 = -pkin(7) * qJDD(4) + (t128 - t167) * qJD(4);
t105 = -qJDD(4) * t75 + (t76 * t92 - t130) * qJD(4);
t81 = t100 * qJDD(1);
t6 = qJDD(4) * pkin(4) - qJD(4) * t24 - t126 * t97 + t81;
t104 = t31 * t28 + t24 * t143 + g(2) * t162 + g(1) * t160 + g(3) * t84 + (-t24 * t91 - t6) * t96;
t7 = qJD(4) * t23 + t97 * qJDD(1) + t100 * t126;
t103 = g(1) * t161 + g(2) * t163 - g(3) * t85 + qJD(5) * t118 - t31 * t30 + t6 * t99 - t96 * t7;
t87 = t92 ^ 2;
t80 = cos(t90);
t79 = sin(t90);
t58 = t77 - t165;
t55 = qJDD(4) * t100 - t102 * t97;
t54 = qJDD(4) * t97 + t100 * t102;
t49 = pkin(2) * t144 + t133;
t48 = t100 * t123;
t47 = t97 * t123;
t34 = 0.2e1 * t127 * t97 + t89 * t93;
t26 = t100 * t122 - t130 * t97;
t25 = t100 * t130 + t122 * t97;
t22 = -0.2e1 * qJD(4) * t151 * t92 + 0.2e1 * t148 * t89;
t12 = t113 * t88 - t21 * t91;
t10 = -t28 ^ 2 + t30 ^ 2;
t9 = t21 * t92 + t120;
t3 = t28 * t91 + t8;
t2 = t20 * t30 + t40 * t8;
t1 = t113 * t8 - t20 * t28 - t21 * t30 - t40 * t9;
t4 = [t140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t54, 0, 0, 0, 0, 0, t12, -t11; 0, qJDD(2), g(1) * t79 - g(2) * t80, g(1) * t80 + g(2) * t79, t89, (t101 * t89 - t129) * pkin(2) + t112, ((-qJDD(2) - t89) * t98 + (-qJD(2) - t92) * t138) * pkin(2) + t153, t34, t22, t54, t55, 0, t105 * t97 + (-t107 - t171) * t100 + t154, t105 * t100 + (t107 - t68) * t97 + t134, t2, t1, t11, t12, 0, t49 * t28 + t58 * t9 + (-qJD(5) * t116 - t96 * t25 + t99 * t26) * t91 + t117 * t88 + t110, t49 * t30 + t58 * t8 - (qJD(5) * t117 + t99 * t25 + t96 * t26) * t91 - t116 * t88 + t111; 0, 0, 0, 0, t89, t112 + t121, (-t137 + (-qJD(3) + t92) * t139) * pkin(2) + t153, t34, t22, t54, t55, 0, t106 * t97 + (-t108 - t171) * t100 + t154, t106 * t100 + (t108 - t68) * t97 + t134, t2, t1, t11, t12, 0, t28 * t133 + t77 * t9 + (-qJD(5) * t114 + t96 * t47 - t99 * t48) * t91 + t115 * t88 + (t147 * t40 - t98 * t28) * t150 + t110, t30 * t133 + t77 * t8 - (qJD(5) * t115 - t99 * t47 - t96 * t48) * t91 - t114 * t88 + (t113 * t147 - t98 * t30) * t150 + t111; 0, 0, 0, 0, 0, 0, 0, -t87 * t148, t151 * t87, t157, t149, qJDD(4), -g(3) * t100 + t109 * t97 + t81, t100 * t109 - t140 * t97, t164, t10, t3, -t120, t88, -(-t23 * t96 - t156) * t91 + (-t143 * t91 - t159 * t28 + t88 * t99) * pkin(4) + t103, (-qJD(5) * t19 + t23 * t91 - t7) * t99 + (-qJD(5) * t91 * t99 - t159 * t30 - t88 * t96) * pkin(4) + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t10, t3, -t120, t88, -t118 * t91 + t103, (-t7 + (-qJD(5) + t91) * t19) * t99 + t104;];
tau_reg = t4;
