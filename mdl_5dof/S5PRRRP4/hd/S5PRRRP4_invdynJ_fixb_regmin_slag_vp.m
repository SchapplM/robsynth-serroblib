% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP4
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
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:23
% EndTime: 2019-12-05 16:46:27
% DurationCPUTime: 1.30s
% Computational Cost: add. (1370->211), mult. (2278->252), div. (0->0), fcn. (1562->10), ass. (0->149)
t88 = qJ(2) + qJ(3);
t82 = cos(t88);
t182 = g(3) * t82;
t89 = sin(pkin(8));
t90 = cos(pkin(8));
t120 = g(1) * t90 + g(2) * t89;
t81 = sin(t88);
t195 = t120 * t81;
t197 = -t182 + t195;
t196 = 2 * qJD(4);
t76 = g(3) * t81;
t135 = t120 * t82 + t76;
t91 = sin(qJ(4));
t94 = cos(qJ(4));
t113 = t94 * pkin(4) + t91 * qJ(5) + pkin(3);
t84 = qJDD(2) + qJDD(3);
t194 = t113 * t84;
t85 = qJD(2) + qJD(3);
t193 = t113 * t85;
t139 = qJD(1) * qJD(2);
t93 = sin(qJ(2));
t96 = cos(qJ(2));
t108 = -t93 * qJDD(1) - t96 * t139;
t138 = qJDD(4) * qJ(5);
t153 = qJD(1) * t93;
t144 = t96 * qJD(1);
t68 = qJD(2) * pkin(2) + t144;
t92 = sin(qJ(3));
t95 = cos(qJ(3));
t32 = t95 * t153 + t92 * t68;
t26 = t85 * pkin(7) + t32;
t164 = t91 * t26;
t80 = t96 * qJDD(1);
t43 = qJDD(2) * pkin(2) - t93 * t139 + t80;
t186 = -(qJD(3) * t68 - t108) * t95 - t92 * t43;
t128 = qJD(3) * t153;
t63 = t92 * t128;
t11 = t84 * pkin(7) - t186 - t63;
t9 = t94 * t11;
t6 = t138 + t9 + (qJD(5) - t164) * qJD(4);
t147 = qJDD(4) * pkin(4);
t148 = qJD(4) * t94;
t187 = t26 * t148 - t147;
t8 = t91 * t11;
t7 = qJDD(5) + t8 + t187;
t191 = t6 * t94 + t7 * t91;
t46 = t92 * t93 - t95 * t96;
t188 = t46 * t85;
t97 = qJD(4) ^ 2;
t185 = pkin(7) * t97;
t181 = g(3) * t91;
t180 = t84 * pkin(3);
t179 = t85 * pkin(3);
t177 = t95 * pkin(2);
t176 = t188 * t85;
t70 = t92 * t153;
t31 = t95 * t68 - t70;
t175 = t31 * t85;
t174 = t32 * t85;
t77 = t92 * pkin(2) + pkin(7);
t173 = t77 * t97;
t170 = t85 * t91;
t169 = t85 * t94;
t168 = t89 * t91;
t167 = t89 * t94;
t166 = t90 * t91;
t165 = t90 * t94;
t162 = t94 * t26;
t161 = t94 * t84;
t152 = qJD(3) * t92;
t133 = pkin(2) * t152;
t142 = qJD(4) * qJ(5);
t145 = t91 * qJD(5);
t149 = qJD(4) * t91;
t39 = pkin(4) * t149 - t94 * t142 - t145;
t28 = t39 + t133;
t47 = t92 * t96 + t95 * t93;
t41 = t47 * qJD(1);
t160 = t28 - t41;
t159 = t39 - t32;
t158 = t195 * t94;
t86 = t91 ^ 2;
t87 = t94 ^ 2;
t157 = t86 - t87;
t156 = t86 + t87;
t154 = pkin(7) * qJDD(4);
t151 = qJD(3) * t95;
t150 = qJD(4) * t85;
t17 = t142 + t162;
t146 = t17 * qJD(4);
t143 = qJDD(1) - g(3);
t141 = qJDD(4) * t77;
t83 = t85 ^ 2;
t137 = t91 * t83 * t94;
t136 = -t82 * t181 + t195 * t91;
t134 = pkin(2) * t151;
t132 = t85 * t152;
t122 = -t68 * t152 + (-t128 + t43) * t95 + t108 * t92;
t12 = -t122 - t180;
t131 = -t12 - t182;
t130 = t156 * t84;
t127 = t31 * t149 + t32 * t169 + t158;
t42 = t95 * t144 - t70;
t126 = t42 * t149 + t41 * t169 + t158;
t125 = t63 + t135;
t124 = -qJD(4) * pkin(4) + qJD(5);
t25 = -t31 - t179;
t123 = t12 * t91 + t25 * t148 - t136;
t121 = -t180 + t185;
t119 = pkin(4) * t91 - qJ(5) * t94;
t16 = t124 + t164;
t118 = t16 * t91 + t17 * t94;
t19 = t85 * t47;
t117 = -t19 * t85 - t46 * t84;
t78 = -pkin(3) - t177;
t116 = t78 * t84 + t173;
t36 = t82 * t167 - t166;
t38 = t82 * t165 + t168;
t115 = g(1) * t38 + g(2) * t36 - t9;
t4 = (t119 * qJD(4) - t145) * t85 - t194 - t122;
t114 = -t185 - t4 + t194;
t44 = -t113 - t177;
t112 = -t44 * t84 - t173 - t4;
t111 = t44 * t85 - t134;
t110 = t78 * t85 - t134;
t35 = t82 * t168 + t165;
t37 = t82 * t166 - t167;
t109 = g(1) * t37 + g(2) * t35 + t81 * t181 - t8;
t107 = t47 * t97 - t117;
t106 = -qJDD(5) + t109;
t105 = -qJDD(4) * t47 + t188 * t196;
t104 = t122 + t197;
t103 = -t91 * t146 + t16 * t148 - t135 + t191;
t102 = -g(3) * t96 + t120 * t93;
t100 = (t16 * t94 - t17 * t91) * qJD(4) + t191;
t99 = -t135 * pkin(7) + t197 * t113;
t98 = qJD(2) ^ 2;
t71 = t91 * t84;
t54 = qJDD(4) * t94 - t97 * t91;
t53 = qJDD(4) * t91 + t97 * t94;
t40 = t119 * t85;
t33 = 0.2e1 * t148 * t170 + t86 * t84;
t23 = -0.2e1 * t157 * t150 + 0.2e1 * t91 * t161;
t21 = t25 * t149;
t14 = -t31 - t193;
t13 = t14 * t149;
t2 = t105 * t91 - t107 * t94;
t1 = t105 * t94 + t107 * t91;
t3 = [t143, 0, t96 * qJDD(2) - t98 * t93, -qJDD(2) * t93 - t98 * t96, 0, t117, -t47 * t84 + t176, 0, 0, 0, 0, 0, t2, t1, t2, t47 * t130 - t156 * t176, -t1, t100 * t47 - t118 * t188 + t14 * t19 + t4 * t46 - g(3); 0, qJDD(2), t80 + t102, t120 * t96 - t143 * t93, t84, t41 * t85 + (t84 * t95 - t132) * pkin(2) + t104, t42 * t85 + (-pkin(2) * t84 - t43) * t92 + ((-pkin(2) * t85 - t68) * qJD(3) + t108) * t95 + t125, t33, t23, t53, t54, 0, t21 + (t110 * qJD(4) - t141) * t91 + (-pkin(2) * t132 - t116 + t131) * t94 + t126, (-t141 + (t110 + t42) * qJD(4)) * t94 + ((-t41 + t133) * t85 + t116) * t91 + t123, t13 + (t111 * qJD(4) - t141) * t91 + (-t28 * t85 + t112 - t182) * t94 + t126, t77 * t130 + t103 + (t134 - t42) * t85 * t156, (t141 + (-t111 - t14 - t42) * qJD(4)) * t94 + (-t160 * t85 + t112) * t91 + t136, t4 * t44 - t118 * t42 + t160 * t14 + (t118 * t151 + t102) * pkin(2) + t100 * t77 + t99; 0, 0, 0, 0, t84, t104 + t174, t125 + t175 + t186, t33, t23, t53, t54, 0, t21 + (-pkin(3) * t150 - t154) * t91 + (-t121 + t131) * t94 + t127, (-t154 + (t31 - t179) * qJD(4)) * t94 + (t121 - t174) * t91 + t123, t13 + (-t113 * t150 - t154) * t91 + (-t39 * t85 + t114 - t182) * t94 + t127, pkin(7) * t130 - t156 * t175 + t103, (t154 + (-t14 - t31 + t193) * qJD(4)) * t94 + (-t159 * t85 + t114) * t91 + t136, t100 * pkin(7) - t113 * t4 - t118 * t31 + t159 * t14 + t99; 0, 0, 0, 0, 0, 0, 0, -t137, t157 * t83, t71, t161, qJDD(4), -t25 * t170 + t109, (-t25 * t85 + t76) * t94 + t115, 0.2e1 * t147 + (-t14 * t91 + t40 * t94) * t85 + t106, -t119 * t84 + ((t17 - t142) * t91 + (t124 - t16) * t94) * t85, -t94 * t76 + 0.2e1 * t138 + qJD(5) * t196 + (t14 * t94 + t40 * t91) * t85 - t115, t6 * qJ(5) - t7 * pkin(4) - t14 * t40 - t16 * t162 - g(1) * (-t37 * pkin(4) + t38 * qJ(5)) - g(2) * (-t35 * pkin(4) + t36 * qJ(5)) + t119 * t76 + (qJD(5) + t164) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t137, t71, -t86 * t83 - t97, t14 * t170 - t106 - t146 + t187;];
tau_reg = t3;
