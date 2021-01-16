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
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 15:42:09
% EndTime: 2021-01-15 15:42:15
% DurationCPUTime: 1.20s
% Computational Cost: add. (1400->218), mult. (3108->300), div. (0->0), fcn. (2285->12), ass. (0->126)
t116 = sin(qJ(5));
t118 = cos(qJ(5));
t113 = sin(pkin(9));
t114 = cos(pkin(9));
t119 = cos(qJ(3));
t154 = t114 * t119;
t145 = qJD(2) * t154;
t117 = sin(qJ(3));
t151 = qJD(2) * t117;
t70 = t113 * t151 - t145;
t77 = t113 * t119 + t114 * t117;
t73 = t77 * qJD(2);
t131 = t116 * t70 - t118 * t73;
t147 = t119 * qJDD(2);
t148 = t117 * qJDD(2);
t136 = t113 * t148 - t114 * t147;
t72 = t77 * qJD(3);
t40 = qJD(2) * t72 + t136;
t149 = qJD(2) * qJD(3);
t144 = t117 * t149;
t125 = t77 * qJDD(2) - t113 * t144;
t143 = t119 * t149;
t41 = t114 * t143 + t125;
t122 = t131 * qJD(5) - t116 * t41 - t118 * t40;
t109 = qJD(3) + qJD(5);
t156 = t131 * t109;
t174 = t122 - t156;
t61 = t118 * t70;
t35 = -t116 * t73 - t61;
t155 = t35 * t109;
t150 = qJD(5) * t116;
t5 = -qJD(5) * t61 - t116 * t40 + t118 * t41 - t73 * t150;
t173 = t5 - t155;
t172 = t131 * t35;
t108 = pkin(8) + qJ(2);
t100 = sin(t108);
t102 = cos(t108);
t138 = g(1) * t102 + g(2) * t100;
t171 = t131 ^ 2 - t35 ^ 2;
t166 = t70 * pkin(7);
t159 = qJ(4) + pkin(6);
t89 = t159 * t119;
t67 = t117 * qJD(1) + qJD(2) * t89;
t157 = t114 * t67;
t158 = qJD(3) * pkin(3);
t88 = t159 * t117;
t65 = t119 * qJD(1) - qJD(2) * t88;
t60 = t65 + t158;
t24 = t113 * t60 + t157;
t13 = t24 - t166;
t99 = t119 * pkin(3) + pkin(2);
t85 = -t99 * qJD(2) + qJD(4);
t44 = t70 * pkin(4) + t85;
t110 = qJ(3) + pkin(9);
t105 = qJ(5) + t110;
t96 = sin(t105);
t97 = cos(t105);
t170 = g(3) * t96 + t13 * t150 + t138 * t97 - t44 * t35;
t104 = t119 * qJDD(1);
t141 = t159 * qJD(3);
t139 = qJD(2) * t141;
t167 = qJD(3) * qJD(1) + qJD(2) * qJD(4) + t159 * qJDD(2);
t22 = qJDD(3) * pkin(3) - t167 * t117 - t119 * t139 + t104;
t25 = (qJDD(1) - t139) * t117 + t167 * t119;
t7 = -t113 * t25 + t114 * t22;
t2 = qJDD(3) * pkin(4) - t41 * pkin(7) + t7;
t8 = t113 * t22 + t114 * t25;
t3 = -t40 * pkin(7) + t8;
t169 = -g(3) * t97 - t116 * t3 + t118 * t2 + t44 * t131 + t138 * t96;
t168 = qJD(5) - t109;
t165 = t73 * pkin(7);
t164 = pkin(3) * t113;
t163 = pkin(3) * t117;
t160 = g(3) * t119;
t53 = t113 * t67;
t28 = t114 * t65 - t53;
t66 = t119 * qJD(4) - t117 * t141;
t68 = -t117 * qJD(4) - t119 * t141;
t29 = t113 * t68 + t114 * t66;
t46 = -t113 * t88 + t114 * t89;
t153 = qJDD(1) - g(3);
t111 = t117 ^ 2;
t152 = -t119 ^ 2 + t111;
t146 = t117 * t158;
t23 = t114 * t60 - t53;
t26 = -t113 * t65 - t157;
t27 = -t113 * t66 + t114 * t68;
t45 = -t113 * t89 - t114 * t88;
t137 = g(1) * t100 - g(2) * t102;
t107 = qJDD(3) + qJDD(5);
t76 = t113 * t117 - t154;
t43 = -t116 * t76 + t118 * t77;
t42 = t116 * t77 + t118 * t76;
t75 = t76 * qJD(3);
t9 = -t42 * qJD(5) - t116 * t72 - t118 * t75;
t135 = t43 * t107 + t9 * t109;
t12 = qJD(3) * pkin(4) - t165 + t23;
t134 = -t116 * t12 - t118 * t13;
t30 = -t77 * pkin(7) + t45;
t31 = -t76 * pkin(7) + t46;
t133 = -t116 * t31 + t118 * t30;
t132 = t116 * t30 + t118 * t31;
t98 = t114 * pkin(3) + pkin(4);
t129 = t116 * t98 + t118 * t164;
t128 = -t116 * t164 + t118 * t98;
t127 = -0.2e1 * pkin(2) * t149 - pkin(6) * qJDD(3);
t64 = pkin(3) * t144 - t99 * qJDD(2) + qJDD(4);
t120 = qJD(3) ^ 2;
t124 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t120 + t137;
t121 = qJD(2) ^ 2;
t123 = t121 * pkin(2) - qJDD(2) * pkin(6) + t138;
t103 = cos(t110);
t101 = sin(t110);
t87 = qJDD(3) * t119 - t120 * t117;
t86 = qJDD(3) * t117 + t120 * t119;
t52 = t76 * pkin(4) - t99;
t48 = t72 * pkin(4) + t146;
t47 = pkin(3) * t151 + t73 * pkin(4);
t18 = t40 * pkin(4) + t64;
t17 = -t72 * pkin(7) + t29;
t16 = t28 - t165;
t15 = t75 * pkin(7) + t27;
t14 = t26 + t166;
t10 = t43 * qJD(5) - t116 * t75 + t118 * t72;
t4 = -t10 * t109 - t42 * t107;
t1 = [t153, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t86, -t72 * qJD(3) - t76 * qJDD(3), t75 * qJD(3) - t77 * qJDD(3), -t77 * t40 + t76 * t41 + t75 * t70 + t72 * t73, -t23 * t72 - t24 * t75 - t7 * t76 + t8 * t77 - g(3), 0, 0, 0, 0, 0, t4, -t135; 0, qJDD(2), t137, t138, t111 * qJDD(2) + 0.2e1 * t117 * t143, 0.2e1 * t117 * t147 - 0.2e1 * t152 * t149, t86, t87, 0, t127 * t117 + t124 * t119, -t124 * t117 + t127 * t119, t45 * qJDD(3) - t99 * t40 + t64 * t76 + t85 * t72 + t137 * t103 + (t70 * t163 + t27) * qJD(3), -t46 * qJDD(3) - t99 * t41 + t64 * t77 - t85 * t75 - t137 * t101 + (t73 * t163 - t29) * qJD(3), t23 * t75 - t24 * t72 - t27 * t73 - t29 * t70 - t46 * t40 - t45 * t41 - t7 * t77 - t8 * t76 - t138, t8 * t46 + t24 * t29 + t7 * t45 + t23 * t27 - t64 * t99 + t85 * t146 - g(1) * (-t100 * t99 + t102 * t159) - g(2) * (t100 * t159 + t102 * t99), -t131 * t9 + t5 * t43, t10 * t131 + t122 * t43 + t35 * t9 - t5 * t42, t135, t4, 0, -t48 * t35 - t52 * t122 + t18 * t42 + t44 * t10 + (-t132 * qJD(5) - t116 * t17 + t118 * t15) * t109 + t133 * t107 + t137 * t97, -t48 * t131 + t52 * t5 + t18 * t43 + t44 * t9 - (t133 * qJD(5) + t116 * t15 + t118 * t17) * t109 - t132 * t107 - t137 * t96; 0, 0, 0, 0, -t117 * t121 * t119, t152 * t121, t148, t147, qJDD(3), t123 * t117 + t104 - t160, -t153 * t117 + t123 * t119, -g(3) * t103 - t26 * qJD(3) - t85 * t73 + t138 * t101 + (qJDD(3) * t114 - t70 * t151) * pkin(3) + t7, g(3) * t101 + t28 * qJD(3) + t85 * t70 + t138 * t103 + (-qJDD(3) * t113 - t73 * t151) * pkin(3) - t8, (t24 + t26) * t73 + (-t23 + t28) * t70 + (-t113 * t40 - t114 * t41) * pkin(3), -t23 * t26 - t24 * t28 + (-t160 + t113 * t8 + t114 * t7 + (-qJD(2) * t85 + t138) * t117) * pkin(3), t172, t171, t173, t174, t107, t128 * t107 + t47 * t35 - (-t116 * t16 + t118 * t14) * t109 + (-t129 * t109 + t134) * qJD(5) + t169, -t129 * t107 - t118 * t3 - t116 * t2 + t47 * t131 + (t116 * t14 + t118 * t16) * t109 + (-t128 * t109 - t118 * t12) * qJD(5) + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t73 * qJD(3) + t136, (-t70 + t145) * qJD(3) + t125, -t70 ^ 2 - t73 ^ 2, t23 * t73 + t24 * t70 - t137 + t64, 0, 0, 0, 0, 0, -t122 - t156, t5 + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t171, t173, t174, t107, t168 * t134 + t169, (-t13 * t109 - t2) * t116 + (-t168 * t12 - t3) * t118 + t170;];
tau_reg = t1;
