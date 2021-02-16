% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:08:52
% EndTime: 2021-01-15 20:08:57
% DurationCPUTime: 0.97s
% Computational Cost: add. (1382->232), mult. (2201->276), div. (0->0), fcn. (1262->12), ass. (0->154)
t180 = qJDD(3) - g(1);
t104 = cos(qJ(2));
t160 = pkin(1) * qJD(2);
t140 = qJD(1) * t160;
t101 = sin(qJ(2));
t144 = qJDD(1) * t101;
t179 = pkin(1) * t144 + t104 * t140;
t96 = qJ(1) + qJ(2);
t82 = pkin(8) + t96;
t71 = sin(t82);
t176 = g(2) * t71;
t72 = cos(t82);
t67 = g(3) * t72;
t178 = t67 - t176;
t103 = cos(qJ(4));
t161 = pkin(1) * qJD(1);
t141 = t104 * t161;
t93 = qJD(1) + qJD(2);
t50 = t93 * pkin(2) + t141;
t142 = t101 * t161;
t98 = cos(pkin(8));
t66 = t98 * t142;
t97 = sin(pkin(8));
t29 = t97 * t50 + t66;
t24 = t93 * pkin(7) + t29;
t133 = qJ(5) * t93 + t24;
t119 = t133 * t103;
t100 = sin(qJ(4));
t94 = t100 ^ 2;
t177 = pkin(4) * t94;
t175 = g(2) * t72;
t92 = qJDD(1) + qJDD(2);
t174 = t92 * pkin(3);
t173 = t98 * pkin(2);
t84 = t103 * qJD(3);
t12 = -t133 * t100 + t84;
t159 = qJD(4) * pkin(4);
t9 = t12 + t159;
t172 = -t12 + t9;
t171 = g(1) * t103;
t170 = t103 * pkin(4);
t169 = t104 * pkin(1);
t41 = t97 * t141 + t66;
t168 = t41 * t93;
t78 = pkin(3) + t170;
t51 = -t78 - t173;
t167 = t51 * t92;
t147 = t100 * qJD(4);
t153 = t103 * t93;
t65 = t97 * t142;
t43 = t98 * t141 - t65;
t166 = t43 * t147 + t41 * t153;
t158 = t100 * t71;
t165 = g(3) * t158 + t100 * t175;
t155 = t101 * t98;
t79 = pkin(2) + t169;
t164 = pkin(1) * t155 + t97 * t79;
t95 = t103 ^ 2;
t163 = -t94 - t95;
t162 = t94 - t95;
t70 = t100 * t92;
t157 = t100 * t93;
t156 = t101 * t97;
t154 = t103 * t92;
t40 = pkin(7) + t164;
t152 = -qJ(5) - t40;
t73 = t97 * pkin(2) + pkin(7);
t151 = -qJ(5) - t73;
t150 = qJD(4) * t93;
t149 = qJDD(4) * pkin(4);
t148 = t100 * t103;
t146 = t103 * qJD(4);
t145 = qJD(4) * qJD(3);
t85 = sin(t96);
t76 = pkin(2) * t85;
t99 = -qJ(5) - pkin(7);
t143 = t71 * t78 + t72 * t99 + t76;
t80 = qJDD(1) * t169;
t37 = t92 * pkin(2) - t101 * t140 + t80;
t17 = t179 * t98 + t97 * t37;
t139 = t93 * t147;
t86 = cos(t96);
t137 = g(2) * t85 - g(3) * t86;
t28 = t98 * t50 - t65;
t136 = -pkin(1) * t156 + t98 * t79;
t18 = -t78 * t93 + qJD(5) - t28;
t16 = -t179 * t97 + t98 * t37;
t62 = pkin(4) * t139;
t111 = qJDD(5) - t16 + t62;
t5 = -t78 * t92 + t111;
t135 = t5 * t100 + t18 * t146 + t165;
t10 = -t16 - t174;
t23 = -t93 * pkin(3) - t28;
t134 = t10 * t100 + t23 * t146 + t165;
t132 = qJD(4) * t152;
t131 = qJD(4) * t151;
t130 = 0.2e1 * t93 * t146;
t129 = qJD(1) * (-qJD(2) + t93);
t128 = qJD(2) * (-qJD(1) - t93);
t11 = t92 * pkin(7) + t17;
t112 = -qJ(5) * t92 - t11 - t145;
t108 = qJD(5) * t93 - t112;
t21 = t24 * t147;
t3 = -t21 + (-qJ(5) * t150 + qJDD(3)) * t100 + t108 * t103;
t127 = t3 * t103 + t178;
t39 = -pkin(3) - t136;
t81 = t103 * qJDD(3);
t125 = g(2) * t158 - t171 + t81;
t77 = pkin(2) * t86;
t124 = -t71 * t99 + t72 * t78 + t77;
t123 = g(3) * t71 + t175;
t122 = -g(2) * t86 - g(3) * t85;
t42 = (t104 * t97 + t155) * t160;
t30 = pkin(4) * t147 + t42;
t34 = t39 - t170;
t121 = t30 * t93 + t34 * t92;
t106 = qJD(4) ^ 2;
t74 = -pkin(3) - t173;
t120 = t106 * t73 + t74 * t92;
t118 = -t123 - t5;
t117 = -t10 - t123;
t116 = t122 + t80;
t115 = -t23 * t93 - t11 - t67;
t114 = -t180 * t100 + t103 * t176 + t21;
t113 = -qJDD(4) * t73 + t74 * t150;
t110 = t106 * t40 + t39 * t92 + t42 * t93;
t44 = (t104 * t98 - t156) * t160;
t109 = -qJDD(4) * t40 + (t39 * t93 - t44) * qJD(4);
t107 = -t67 + (-qJD(5) - t18) * t93 + t112;
t105 = cos(qJ(1));
t102 = sin(qJ(1));
t91 = t93 ^ 2;
t89 = t105 * pkin(1);
t88 = t102 * pkin(1);
t87 = t103 * qJ(5);
t83 = t103 * qJD(5);
t53 = qJDD(4) * t103 - t106 * t100;
t52 = qJDD(4) * t100 + t106 * t103;
t47 = t103 * t73 + t87;
t46 = t151 * t100;
t38 = t100 * t130 + t94 * t92;
t36 = -t100 * qJD(5) + t103 * t131;
t35 = t100 * t131 + t83;
t32 = t43 * t146;
t26 = t103 * t40 + t87;
t25 = t152 * t100;
t22 = 0.2e1 * t92 * t148 - 0.2e1 * t162 * t150;
t19 = t23 * t147;
t14 = t18 * t147;
t13 = t100 * qJD(3) + t119;
t7 = (-qJD(5) - t44) * t100 + t103 * t132;
t6 = t100 * t132 + t103 * t44 + t83;
t2 = -qJD(4) * t119 - t108 * t100 + t149 + t81;
t1 = [qJDD(1), -g(2) * t105 - g(3) * t102, g(2) * t102 - g(3) * t105, t92, (t101 * t128 + t104 * t92) * pkin(1) + t116, ((-qJDD(1) - t92) * t101 + t104 * t128) * pkin(1) + t137, t17 * t164 + t29 * t44 + t16 * t136 - t28 * t42 - g(2) * (t77 + t89) - g(3) * (t76 + t88), t38, t22, t52, t53, 0, t19 + t109 * t100 + (-t110 + t117) * t103, t100 * t110 + t103 * t109 + t134, t25 * qJDD(4) + t14 + (t34 * t157 + t7) * qJD(4) + (t118 - t121) * t103, -t26 * qJDD(4) + t121 * t100 + (t34 * t153 - t6) * qJD(4) + t135, (t26 * t92 + t6 * t93 + (-t25 * t93 - t9) * qJD(4)) * t103 + (-t25 * t92 - t7 * t93 - t2 + (-t26 * t93 - t13) * qJD(4)) * t100 + t127, t3 * t26 + t13 * t6 + t2 * t25 + t9 * t7 + t5 * t34 + t18 * t30 - g(2) * (t124 + t89) - g(3) * (t88 + t143); 0, 0, 0, t92, t101 * pkin(1) * t129 + t116, (t104 * t129 - t144) * pkin(1) + t137, t28 * t41 - t29 * t43 + (t16 * t98 + t17 * t97 + t122) * pkin(2), t38, t22, t52, t53, 0, t19 + t113 * t100 + (t117 - t120) * t103 + t166, t32 + t113 * t103 + (t120 - t168) * t100 + t134, t46 * qJDD(4) + t14 + (t51 * t157 + t36) * qJD(4) + (t118 - t62 - t167) * t103 + t166, -t47 * qJDD(4) + t32 + (t167 - t168) * t100 + (-t35 + (t103 * t51 + t177) * t93) * qJD(4) + t135, (-qJD(4) * t9 + t47 * t92) * t103 + (-t13 * qJD(4) - t46 * t92 - t2) * t100 + (-t100 * t36 + t103 * t35 + t163 * t43 + (-t100 * t47 - t103 * t46) * qJD(4)) * t93 + t127, t3 * t47 + t2 * t46 + t9 * t36 + t5 * t51 - t18 * t41 - g(2) * t124 - g(3) * t143 + (-t103 * t43 + t35) * t13 + (t18 * t159 + t9 * t43) * t100; 0, 0, 0, 0, 0, 0, t180, 0, 0, 0, 0, 0, t53, -t52, t53, -t52, 0, t3 * t100 + t2 * t103 - g(1) + (-t9 * t100 + t13 * t103) * qJD(4); 0, 0, 0, 0, 0, 0, 0, -t91 * t148, t162 * t91, t70, t154, qJDD(4), t100 * t115 + t125, (-t100 * t24 + t84) * qJD(4) + (t115 - t145) * t103 + t114, 0.2e1 * t149 + (t13 - t119) * qJD(4) + (t91 * t170 + t107) * t100 + t125, -t91 * t177 + (qJ(5) * t157 + t12) * qJD(4) + t107 * t103 + t114, -pkin(4) * t70 + (-t159 + t172) * t153, t172 * t13 + (-t171 + t2 + (-t18 * t93 - t178) * t100) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t139 - t154, t70 + t130, t163 * t91, t9 * t157 - t174 + (-pkin(4) * t92 - t13 * t93) * t103 + t111 + t123;];
tau_reg = t1;
