% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:37
% EndTime: 2019-12-05 18:30:39
% DurationCPUTime: 1.19s
% Computational Cost: add. (2923->218), mult. (4964->265), div. (0->0), fcn. (2944->16), ass. (0->148)
t98 = qJ(1) + qJ(2);
t90 = pkin(9) + t98;
t83 = qJ(4) + t90;
t74 = sin(t83);
t75 = cos(t83);
t135 = g(2) * t75 + g(3) * t74;
t102 = sin(qJ(4));
t106 = cos(qJ(4));
t107 = cos(qJ(2));
t100 = cos(pkin(9));
t103 = sin(qJ(2));
t154 = t100 * t103;
t99 = sin(pkin(9));
t127 = pkin(1) * (-t107 * t99 - t154);
t50 = qJD(1) * t127;
t159 = t103 * t99;
t128 = pkin(1) * (t100 * t107 - t159);
t52 = qJD(1) * t128;
t181 = pkin(2) * t99;
t170 = pkin(2) * t100;
t79 = pkin(3) + t170;
t55 = -t102 * t181 + t106 * t79;
t165 = t55 * qJD(4) - t102 * t50 - t106 * t52;
t94 = qJD(1) + qJD(2);
t88 = qJD(4) + t94;
t188 = t165 * t88;
t150 = qJDD(1) * t103;
t152 = qJD(1) * t107;
t121 = pkin(1) * (qJD(2) * t152 + t150);
t187 = t100 * t121;
t109 = qJD(5) ^ 2;
t57 = t102 * t79 + t106 * t181;
t164 = t57 * qJD(4) - t102 * t52 + t106 * t50;
t144 = t164 * t88;
t48 = -pkin(4) - t55;
t49 = pkin(8) + t57;
t93 = qJDD(1) + qJDD(2);
t87 = qJDD(4) + t93;
t186 = t109 * t49 + t48 * t87 + t144;
t174 = pkin(1) * t103;
t147 = qJD(1) * t174;
t61 = pkin(1) * t152 + pkin(2) * t94;
t35 = t100 * t61 - t147 * t99;
t33 = pkin(3) * t94 + t35;
t36 = t100 * t147 + t61 * t99;
t17 = -t102 * t36 + t106 * t33;
t101 = sin(qJ(5));
t105 = cos(qJ(5));
t18 = t102 * t33 + t106 * t36;
t16 = pkin(8) * t88 + t18;
t11 = qJD(3) * t105 - t101 * t16;
t158 = t105 * t16;
t12 = qJD(3) * t101 + t158;
t156 = qJD(5) * t11;
t182 = pkin(2) * t93;
t172 = pkin(1) * t107;
t85 = qJDD(1) * t172;
t46 = -qJD(2) * t147 + t182 + t85;
t29 = t100 * t46 - t121 * t99;
t21 = pkin(3) * t93 + t29;
t166 = t46 * t99;
t30 = t166 + t187;
t141 = -qJD(4) * t17 - t102 * t21 - t106 * t30;
t5 = pkin(8) * t87 - t141;
t2 = qJDD(3) * t101 + t105 * t5 + t156;
t155 = qJD(5) * t12;
t89 = t105 * qJDD(3);
t3 = -t101 * t5 - t155 + t89;
t115 = -t3 * t101 + t2 * t105 + (-t101 * t12 - t105 * t11) * qJD(5);
t69 = g(2) * t74;
t163 = -g(3) * t75 + t69;
t91 = sin(t98);
t92 = cos(t98);
t185 = g(2) * t92 + g(3) * t91;
t8 = -qJD(4) * t18 - t102 * t30 + t106 * t21;
t184 = pkin(2) * t91;
t183 = pkin(2) * t92;
t77 = sin(t90);
t180 = pkin(3) * t77;
t78 = cos(t90);
t179 = pkin(3) * t78;
t178 = pkin(4) * t87;
t177 = pkin(4) * t88;
t84 = pkin(2) + t172;
t54 = -pkin(1) * t159 + t100 * t84;
t47 = pkin(3) + t54;
t56 = pkin(1) * t154 + t84 * t99;
t27 = -t102 * t56 + t106 * t47;
t51 = qJD(2) * t127;
t53 = qJD(2) * t128;
t9 = qJD(4) * t27 + t102 * t51 + t106 * t53;
t176 = t88 * t9;
t15 = -t17 - t177;
t151 = qJD(5) * t105;
t6 = -t178 - t8;
t175 = t6 * t101 + t15 * t151;
t104 = sin(qJ(1));
t173 = pkin(1) * t104;
t108 = cos(qJ(1));
t171 = pkin(1) * t108;
t28 = t102 * t47 + t106 * t56;
t10 = t28 * qJD(4) + t102 * t53 - t106 * t51;
t169 = t10 * t88;
t168 = t17 * t88;
t167 = t18 * t88;
t95 = t101 ^ 2;
t96 = t105 ^ 2;
t162 = t95 - t96;
t161 = t95 + t96;
t153 = t101 * t105;
t97 = qJDD(3) - g(1);
t149 = t15 * qJD(5) * t101 + t105 * t135;
t148 = t85 + t185;
t86 = t88 ^ 2;
t146 = t86 * t153;
t145 = -g(2) * t91 + g(3) * t92;
t143 = t161 * t87;
t140 = qJD(1) * (-qJD(2) + t94);
t139 = qJD(2) * (-qJD(1) - t94);
t138 = t101 * t88 * t151;
t137 = -t180 - t184;
t136 = -t179 - t183;
t134 = -t173 - t184;
t133 = -t171 - t183;
t132 = g(2) * t108 + g(3) * t104;
t130 = t101 * t11 - t105 * t12;
t129 = t141 - t163;
t125 = -pkin(8) * t109 + t167 + t178;
t23 = -pkin(4) - t27;
t24 = pkin(8) + t28;
t124 = -t109 * t24 - t23 * t87 - t169;
t66 = t75 * pkin(8);
t123 = -pkin(4) * t74 + t137 + t66;
t122 = -pkin(8) * qJDD(5) + (t17 - t177) * qJD(5);
t120 = -qJDD(5) * t24 + (t23 * t88 - t9) * qJD(5);
t119 = -pkin(4) * t75 - pkin(8) * t74 + t136;
t118 = -qJDD(5) * t49 + (t48 * t88 - t165) * qJD(5);
t116 = -qJD(3) * qJD(5) - t15 * t88 - t163 - t5;
t114 = t163 + t115;
t113 = g(2) * t78 + g(3) * t77 + t29;
t112 = t8 + t135;
t111 = -g(2) * t77 + g(3) * t78 - t187;
t63 = qJDD(5) * t105 - t101 * t109;
t62 = qJDD(5) * t101 + t105 * t109;
t39 = t87 * t96 - 0.2e1 * t138;
t38 = t87 * t95 + 0.2e1 * t138;
t32 = -0.2e1 * qJD(5) * t162 * t88 + 0.2e1 * t153 * t87;
t1 = [0, 0, 0, 0, 0, qJDD(1), t132, -g(2) * t104 + g(3) * t108, 0, 0, 0, 0, 0, 0, 0, t93, (t103 * t139 + t107 * t93) * pkin(1) + t148, ((-qJDD(1) - t93) * t103 + t107 * t139) * pkin(1) + t145, 0, (t132 + (t103 ^ 2 + t107 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t93, t51 * t94 + t54 * t93 + t113, -t53 * t94 - t56 * t93 + t111 - t166, 0, -g(2) * t133 - g(3) * t134 + t29 * t54 + t30 * t56 + t35 * t51 + t36 * t53, 0, 0, 0, 0, 0, t87, t27 * t87 + t112 - t169, -t28 * t87 + t129 - t176, 0, -t141 * t28 + t18 * t9 + t8 * t27 - t17 * t10 - g(2) * (t133 - t179) - g(3) * (t134 - t180), t38, t32, t62, t39, t63, 0, t120 * t101 + (t124 - t6) * t105 + t149, t120 * t105 + (-t124 - t135) * t101 + t175, t143 * t24 + t161 * t176 + t114, t6 * t23 + t15 * t10 - g(2) * (t119 - t171) - g(3) * (t123 - t173) - t130 * t9 + t115 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t140 * t174 + t148, (t107 * t140 - t150) * pkin(1) + t145, 0, 0, 0, 0, 0, 0, 0, t93, t93 * t170 - t50 * t94 + t113, t52 * t94 + (-t46 - t182) * t99 + t111, 0, -t35 * t50 - t36 * t52 + (t100 * t29 + t30 * t99 + t185) * pkin(2), 0, 0, 0, 0, 0, t87, t55 * t87 + t112 - t144, -t57 * t87 + t129 - t188, 0, -g(2) * t136 - g(3) * t137 - t141 * t57 - t164 * t17 + t165 * t18 + t8 * t55, t38, t32, t62, t39, t63, 0, t118 * t101 + (-t6 - t186) * t105 + t149, t118 * t105 + (-t135 + t186) * t101 + t175, t49 * t143 + t161 * t188 + t114, t6 * t48 - g(2) * t119 - g(3) * t123 + t164 * t15 + ((t2 - t156) * t49 + t165 * t12) * t105 + ((-t3 - t155) * t49 - t165 * t11) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, t63, -t62, 0, -qJD(5) * t130 + t101 * t2 + t105 * t3 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t112 + t167, t129 + t168, 0, 0, t38, t32, t62, t39, t63, 0, t122 * t101 + (t125 - t6) * t105 + t149, t122 * t105 + (-t125 - t135) * t101 + t175, pkin(8) * t143 - t161 * t168 + t114, -g(3) * t66 - t15 * t18 + t130 * t17 + (t135 - t6) * pkin(4) + (t115 + t69) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, t162 * t86, t101 * t87, t146, t105 * t87, qJDD(5), -g(1) * t105 + t89 + (t12 - t158) * qJD(5) + t116 * t101, t156 + (qJD(5) * t16 - t97) * t101 + t116 * t105, 0, 0;];
tau_reg = t1;
