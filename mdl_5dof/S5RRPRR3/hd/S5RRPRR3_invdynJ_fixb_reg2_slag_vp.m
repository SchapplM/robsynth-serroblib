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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:34:20
% EndTime: 2022-01-20 10:34:25
% DurationCPUTime: 1.38s
% Computational Cost: add. (2923->218), mult. (4964->266), div. (0->0), fcn. (2944->16), ass. (0->150)
t102 = qJ(1) + qJ(2);
t93 = pkin(9) + t102;
t86 = qJ(4) + t93;
t77 = cos(t86);
t186 = g(2) * t77;
t97 = qJDD(1) + qJDD(2);
t90 = qJDD(4) + t97;
t183 = t90 * pkin(4);
t106 = sin(qJ(4));
t110 = cos(qJ(4));
t103 = sin(pkin(9));
t104 = cos(pkin(9));
t107 = sin(qJ(2));
t182 = pkin(1) * t107;
t151 = qJD(1) * t182;
t111 = cos(qJ(2));
t155 = qJD(1) * t111;
t98 = qJD(1) + qJD(2);
t61 = pkin(1) * t155 + pkin(2) * t98;
t35 = -t103 * t151 + t104 * t61;
t33 = pkin(3) * t98 + t35;
t36 = t103 * t61 + t104 * t151;
t18 = t106 * t33 + t110 * t36;
t153 = qJDD(1) * t107;
t124 = pkin(1) * (qJD(2) * t155 + t153);
t189 = pkin(2) * t97;
t177 = t111 * pkin(1);
t88 = qJDD(1) * t177;
t46 = -qJD(2) * t151 + t189 + t88;
t29 = -t103 * t124 + t104 * t46;
t21 = t97 * pkin(3) + t29;
t166 = t103 * t46;
t194 = t104 * t124;
t30 = t166 + t194;
t8 = -t18 * qJD(4) - t106 * t30 + t110 * t21;
t6 = -t183 - t8;
t195 = t6 + t186;
t113 = qJD(5) ^ 2;
t159 = t104 * t107;
t131 = pkin(1) * (-t103 * t111 - t159);
t50 = qJD(1) * t131;
t160 = t103 * t107;
t130 = pkin(1) * (t104 * t111 - t160);
t52 = qJD(1) * t130;
t180 = pkin(2) * t103;
t178 = t104 * pkin(2);
t81 = pkin(3) + t178;
t57 = t106 * t81 + t110 * t180;
t172 = t57 * qJD(4) - t106 * t52 + t110 * t50;
t91 = qJD(4) + t98;
t147 = t172 * t91;
t55 = -t106 * t180 + t110 * t81;
t48 = -pkin(4) - t55;
t49 = pkin(8) + t57;
t193 = -t113 * t49 - t48 * t90 - t147;
t17 = -t106 * t36 + t110 * t33;
t105 = sin(qJ(5));
t109 = cos(qJ(5));
t16 = pkin(8) * t91 + t18;
t11 = qJD(3) * t109 - t105 * t16;
t164 = t109 * t16;
t12 = qJD(3) * t105 + t164;
t157 = t11 * qJD(5);
t144 = -t17 * qJD(4) - t106 * t21 - t110 * t30;
t5 = pkin(8) * t90 - t144;
t2 = t105 * qJDD(3) + t109 * t5 + t157;
t156 = t12 * qJD(5);
t92 = t109 * qJDD(3);
t3 = -t105 * t5 - t156 + t92;
t118 = -t3 * t105 + t2 * t109 + (-t105 * t12 - t109 * t11) * qJD(5);
t76 = sin(t86);
t170 = -g(1) * t77 - g(2) * t76;
t173 = t55 * qJD(4) - t106 * t50 - t110 * t52;
t192 = t173 * t91;
t100 = t109 ^ 2;
t99 = t105 ^ 2;
t161 = t100 + t99;
t94 = sin(t102);
t95 = cos(t102);
t191 = g(1) * t94 - g(2) * t95;
t190 = pkin(2) * t94;
t79 = sin(t93);
t188 = pkin(3) * t79;
t187 = pkin(4) * t91;
t70 = g(1) * t76;
t87 = pkin(2) + t177;
t54 = -pkin(1) * t160 + t104 * t87;
t47 = pkin(3) + t54;
t56 = pkin(1) * t159 + t103 * t87;
t27 = -t106 * t56 + t110 * t47;
t51 = qJD(2) * t131;
t53 = qJD(2) * t130;
t9 = t27 * qJD(4) + t106 * t51 + t110 * t53;
t184 = t9 * t91;
t108 = sin(qJ(1));
t181 = pkin(1) * t108;
t28 = t106 * t47 + t110 * t56;
t10 = t28 * qJD(4) + t106 * t53 - t110 * t51;
t179 = t10 * t91;
t176 = t17 * t91;
t175 = t18 * t91;
t15 = -t17 - t187;
t174 = t15 * qJD(5) * t105 + t109 * t70;
t171 = t77 * pkin(4) + t76 * pkin(8);
t80 = cos(t93);
t75 = pkin(3) * t80;
t85 = pkin(2) * t95;
t169 = t75 + t85;
t168 = g(1) * t95 + g(2) * t94;
t112 = cos(qJ(1));
t96 = t112 * pkin(1);
t167 = t85 + t96;
t162 = t100 - t99;
t158 = t105 * t109;
t101 = qJDD(3) - g(3);
t154 = qJD(5) * t109;
t152 = t195 * t105 + t15 * t154;
t89 = t91 ^ 2;
t150 = t89 * t158;
t148 = -pkin(4) * t76 + t77 * pkin(8);
t146 = t161 * t90;
t143 = t169 + t171;
t142 = qJD(1) * (-qJD(2) + t98);
t141 = qJD(2) * (-qJD(1) - t98);
t140 = t88 + t191;
t139 = t105 * t91 * t154;
t138 = -t188 - t190;
t137 = -t181 - t190;
t135 = g(1) * t108 - g(2) * t112;
t133 = t105 * t11 - t109 * t12;
t132 = t144 - t170;
t128 = pkin(8) * t113 - t175 - t183;
t23 = -pkin(4) - t27;
t24 = pkin(8) + t28;
t127 = t113 * t24 + t23 * t90 + t179;
t126 = t138 + t148;
t125 = -pkin(8) * qJDD(5) + (t17 - t187) * qJD(5);
t123 = -qJDD(5) * t24 + (t23 * t91 - t9) * qJD(5);
t122 = -qJDD(5) * t49 + (t48 * t91 - t173) * qJD(5);
t120 = -qJD(3) * qJD(5) - t15 * t91 - t170 - t5;
t119 = g(1) * t80 + g(2) * t79 - t194;
t117 = t170 + t118;
t116 = g(1) * t79 - g(2) * t80 + t29;
t115 = t70 + t8 - t186;
t63 = qJDD(5) * t109 - t105 * t113;
t62 = qJDD(5) * t105 + t109 * t113;
t39 = t100 * t90 - 0.2e1 * t139;
t38 = t90 * t99 + 0.2e1 * t139;
t32 = 0.2e1 * t162 * t91 * qJD(5) + 0.2e1 * t90 * t158;
t1 = [0, 0, 0, 0, 0, qJDD(1), t135, g(1) * t112 + g(2) * t108, 0, 0, 0, 0, 0, 0, 0, t97, (t107 * t141 + t111 * t97) * pkin(1) + t140, ((-qJDD(1) - t97) * t107 + t111 * t141) * pkin(1) + t168, 0, (t135 + (t107 ^ 2 + t111 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t97, t51 * t98 + t54 * t97 + t116, -t53 * t98 - t56 * t97 + t119 - t166, 0, -g(1) * t137 - g(2) * t167 + t29 * t54 + t30 * t56 + t35 * t51 + t36 * t53, 0, 0, 0, 0, 0, t90, t27 * t90 + t115 - t179, -t28 * t90 + t132 - t184, 0, -t144 * t28 + t18 * t9 + t8 * t27 - t17 * t10 - g(1) * (t137 - t188) - g(2) * (t75 + t167), t38, t32, t62, t39, t63, 0, t123 * t105 + (-t127 - t195) * t109 + t174, t123 * t109 + (t127 - t70) * t105 + t152, t117 + t161 * (t24 * t90 + t184), t6 * t23 + t15 * t10 - g(1) * (t126 - t181) - g(2) * (t96 + t143) - t133 * t9 + t118 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t142 * t182 + t140, (t111 * t142 - t153) * pkin(1) + t168, 0, 0, 0, 0, 0, 0, 0, t97, t97 * t178 - t50 * t98 + t116, t52 * t98 + (-t46 - t189) * t103 + t119, 0, -t35 * t50 - t36 * t52 + (t103 * t30 + t104 * t29 + t191) * pkin(2), 0, 0, 0, 0, 0, t90, t55 * t90 + t115 - t147, -t57 * t90 + t132 - t192, 0, -g(1) * t138 - g(2) * t169 - t144 * t57 - t172 * t17 + t173 * t18 + t8 * t55, t38, t32, t62, t39, t63, 0, t122 * t105 + (-t195 + t193) * t109 + t174, t122 * t109 + (-t70 - t193) * t105 + t152, t49 * t146 + t161 * t192 + t117, t6 * t48 - g(1) * t126 - g(2) * t143 + t172 * t15 + ((t2 - t157) * t49 + t173 * t12) * t109 + ((-t3 - t156) * t49 - t173 * t11) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, 0, 0, 0, 0, 0, t63, -t62, 0, -t133 * qJD(5) + t2 * t105 + t3 * t109 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t115 + t175, t132 + t176, 0, 0, t38, t32, t62, t39, t63, 0, t125 * t105 + (-t128 - t195) * t109 + t174, t125 * t109 + (t128 - t70) * t105 + t152, pkin(8) * t146 - t161 * t176 + t117, -t6 * pkin(4) + t118 * pkin(8) - g(1) * t148 - g(2) * t171 + t133 * t17 - t15 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, -t162 * t89, t105 * t90, t150, t109 * t90, qJDD(5), -g(3) * t109 + t92 + (t12 - t164) * qJD(5) + t120 * t105, t157 + (qJD(5) * t16 - t101) * t105 + t120 * t109, 0, 0;];
tau_reg = t1;
