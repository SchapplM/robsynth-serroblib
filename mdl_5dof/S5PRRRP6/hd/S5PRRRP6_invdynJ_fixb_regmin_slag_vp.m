% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP6
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
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:21
% EndTime: 2019-12-05 16:52:26
% DurationCPUTime: 1.55s
% Computational Cost: add. (1759->255), mult. (3713->339), div. (0->0), fcn. (2617->10), ass. (0->146)
t102 = sin(qJ(4));
t105 = cos(qJ(3));
t187 = cos(qJ(4));
t142 = qJDD(2) * t187;
t103 = sin(qJ(3));
t153 = t103 * qJDD(2);
t61 = t102 * t105 + t187 * t103;
t96 = qJD(3) + qJD(4);
t198 = t96 * t61;
t16 = qJD(2) * t198 + t102 * t153 - t105 * t142;
t106 = cos(qJ(2));
t100 = sin(pkin(8));
t101 = cos(pkin(8));
t135 = g(1) * t101 + g(2) * t100;
t128 = t135 * t106;
t104 = sin(qJ(2));
t184 = g(3) * t104;
t197 = t128 + t184;
t183 = g(3) * t106;
t148 = t187 * t105;
t168 = t102 * t103;
t127 = t148 - t168;
t123 = t106 * t127;
t189 = pkin(7) + pkin(6);
t66 = t189 * t103;
t67 = t189 * t105;
t129 = -t102 * t67 - t187 * t66;
t149 = qJD(3) * t189;
t62 = t103 * t149;
t63 = t105 * t149;
t188 = -qJD(1) * t123 + t129 * qJD(4) - t102 * t63 - t187 * t62;
t36 = -t102 * t66 + t187 * t67;
t55 = t61 * qJD(2);
t99 = qJ(3) + qJ(4);
t93 = sin(t99);
t95 = qJDD(3) + qJDD(4);
t196 = -t104 * (qJD(1) * t55 + t135 * t93) + t93 * t183 - t188 * t96 - t36 * t95;
t90 = t95 * qJ(5);
t91 = t96 * qJD(5);
t195 = t90 + t91;
t52 = t127 * t104;
t92 = t95 * pkin(4);
t194 = qJDD(5) - t92;
t108 = qJD(3) ^ 2;
t157 = qJD(1) * qJD(2);
t138 = -qJDD(1) * t106 + t104 * t157;
t191 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t108 + (t135 + t157) * t104 - t138 - t183;
t190 = t55 ^ 2;
t163 = qJD(1) * t104;
t73 = qJD(2) * pkin(6) + t163;
t143 = pkin(7) * qJD(2) + t73;
t50 = t143 * t105;
t175 = t102 * t50;
t49 = t143 * t103;
t41 = qJD(3) * pkin(3) - t49;
t22 = t187 * t41 - t175;
t182 = t22 * t96;
t150 = t187 * t50;
t23 = t102 * t41 + t150;
t181 = t23 * t96;
t137 = qJD(2) * t148;
t161 = qJD(2) * t103;
t147 = t102 * t161;
t53 = -t137 + t147;
t180 = t55 * t53;
t124 = t106 * t61;
t179 = -qJD(1) * t124 + t36 * qJD(4) - t102 * t62 + t187 * t63;
t146 = qJD(4) * t187;
t26 = -t187 * t49 - t175;
t178 = pkin(3) * t146 + qJD(5) - t26;
t97 = t103 ^ 2;
t177 = -t105 ^ 2 + t97;
t176 = qJD(2) * pkin(2);
t174 = t104 * t93;
t94 = cos(t99);
t173 = t104 * t94;
t170 = t100 * t106;
t169 = t101 * t106;
t167 = t103 * t106;
t166 = qJD(5) - t22;
t165 = qJDD(1) - g(3);
t109 = qJD(2) ^ 2;
t164 = t108 + t109;
t162 = qJD(1) * t106;
t160 = qJD(2) * t104;
t159 = qJD(3) * t103;
t158 = qJD(4) * t102;
t156 = qJD(2) * qJD(3);
t155 = qJDD(2) * t106;
t154 = qJDD(3) * t103;
t152 = t105 * qJDD(2);
t151 = pkin(3) * t159;
t89 = pkin(3) * t105 + pkin(2);
t145 = t103 * t156;
t144 = t105 * t156;
t58 = qJDD(2) * pkin(6) + qJDD(1) * t104 + t106 * t157;
t21 = -qJD(3) * t105 * t73 + qJDD(3) * pkin(3) - t103 * t58 + (-t144 - t153) * pkin(7);
t24 = -t73 * t159 + t105 * t58 + (-t145 + t152) * pkin(7);
t141 = t102 * t21 + t41 * t146 - t50 * t158 + t187 * t24;
t140 = t102 * t24 + t50 * t146 + t41 * t158 - t187 * t21;
t139 = -t102 * t152 - t103 * t142 - t96 * t137;
t25 = -t102 * t49 + t150;
t136 = -pkin(3) * t158 + t25;
t28 = pkin(4) * t55 + qJ(5) * t53;
t134 = g(1) * t100 - g(2) * t101;
t132 = t96 * t168;
t130 = pkin(4) * t94 + qJ(5) * t93 + t89;
t51 = t61 * t104;
t8 = qJD(2) * t124 + t52 * t96;
t125 = -t106 * t16 + t53 * t160 - t51 * t95 - t8 * t96;
t59 = -t89 * qJD(2) - t162;
t46 = -t101 * t93 + t94 * t170;
t48 = t100 * t93 + t94 * t169;
t121 = g(1) * t48 + g(2) * t46 + g(3) * t173 - t141;
t45 = t101 * t94 + t93 * t170;
t47 = -t100 * t94 + t93 * t169;
t120 = g(1) * t47 + g(2) * t45 + g(3) * t174 - t140;
t32 = pkin(3) * t145 - t89 * qJDD(2) + t138;
t15 = qJD(2) * t132 + t139;
t7 = qJD(2) * t123 - t104 * t198;
t119 = -t106 * t15 - t55 * t160 + t52 * t95 + t7 * t96;
t74 = -t162 - t176;
t118 = -pkin(6) * qJDD(3) + (t162 + t74 - t176) * qJD(3);
t117 = t129 * t95 + t135 * t173 - t179 * t96 - t94 * t183;
t17 = pkin(4) * t53 - qJ(5) * t55 + t59;
t116 = -t17 * t53 - t121;
t115 = t59 * t53 + t121;
t114 = -t59 * t55 + t120;
t113 = -t74 * qJD(2) + t197 - t58;
t112 = -t17 * t55 + t120 - t194;
t111 = -g(1) * (-t47 * pkin(4) + qJ(5) * t48) - g(2) * (-t45 * pkin(4) + qJ(5) * t46) - g(3) * (-pkin(4) * t174 + qJ(5) * t173);
t88 = -t187 * pkin(3) - pkin(4);
t85 = pkin(3) * t102 + qJ(5);
t30 = -qJD(3) * t148 - t105 * t146 + t132;
t29 = -pkin(4) * t127 - qJ(5) * t61 - t89;
t27 = pkin(3) * t161 + t28;
t18 = -t53 ^ 2 + t190;
t12 = t96 * qJ(5) + t23;
t11 = -t96 * pkin(4) + t166;
t6 = t55 * t96 - t16;
t5 = -t139 + (-t147 + t53) * t96;
t4 = pkin(4) * t198 + qJ(5) * t30 - qJD(5) * t61 + t151;
t3 = pkin(4) * t16 + qJ(5) * t15 - qJD(5) * t55 + t32;
t2 = t140 + t194;
t1 = t141 + t195;
t9 = [t165, 0, -t104 * t109 + t155, -qJDD(2) * t104 - t106 * t109, 0, 0, 0, 0, 0, (-0.2e1 * t145 + t152) * t106 + (-t164 * t105 - t154) * t104, (-qJDD(3) * t104 - 0.2e1 * t106 * t156) * t105 + (t164 * t104 - t155) * t103, 0, 0, 0, 0, 0, t125, -t119, t125, -t15 * t51 - t16 * t52 - t53 * t7 + t55 * t8, t119, t1 * t52 - t106 * t3 + t11 * t8 + t12 * t7 + t17 * t160 + t2 * t51 - g(3); 0, qJDD(2), t135 * t104 + t165 * t106, -t165 * t104 + t128, qJDD(2) * t97 + 0.2e1 * t103 * t144, 0.2e1 * t103 * t152 - 0.2e1 * t177 * t156, t105 * t108 + t154, qJDD(3) * t105 - t103 * t108, 0, t118 * t103 + t105 * t191, -t103 * t191 + t118 * t105, -t15 * t61 - t55 * t30, -t127 * t15 - t16 * t61 - t198 * t55 + t30 * t53, -t30 * t96 + t61 * t95, t127 * t95 - t198 * t96, 0, -t16 * t89 + t198 * t59 - t32 * t127 + (t151 - t163) * t53 + t117, t15 * t89 + t55 * t151 - t30 * t59 + t32 * t61 + t196, t16 * t29 + t17 * t198 - t3 * t127 + (t4 - t163) * t53 + t117, t1 * t127 - t11 * t30 - t12 * t198 + t129 * t15 - t16 * t36 + t179 * t55 - t188 * t53 + t2 * t61 - t197, t15 * t29 + t17 * t30 - t3 * t61 - t4 * t55 - t196, t1 * t36 + t17 * t4 - t2 * t129 + t3 * t29 + t188 * t12 + t179 * t11 + (-g(3) * t130 - t135 * t189) * t106 + (-g(3) * t189 - t17 * qJD(1) + t130 * t135) * t104; 0, 0, 0, 0, -t103 * t109 * t105, t177 * t109, t153, t152, qJDD(3), t113 * t103 - t134 * t105, t134 * t103 + t113 * t105, t180, t18, t5, t6, t95, t25 * t96 + (-t96 * t158 - t53 * t161 + t187 * t95) * pkin(3) + t114, t26 * t96 + (-t102 * t95 - t96 * t146 - t55 * t161) * pkin(3) + t115, t136 * t96 - t27 * t53 - t88 * t95 + t112, -t15 * t88 - t16 * t85 + (t12 - t136) * t55 + (t11 - t178) * t53, t178 * t96 + t27 * t55 + t85 * t95 + t116 + t195, t1 * t85 + t2 * t88 - t17 * t27 - t11 * t25 + t178 * t12 + (t11 * t158 - g(1) * (t100 * t105 - t101 * t167) - g(2) * (-t100 * t167 - t101 * t105) + t103 * t184) * pkin(3) + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t18, t5, t6, t95, t114 + t181, t115 + t182, -t28 * t53 + t112 + t181 + t92, pkin(4) * t15 - qJ(5) * t16 + (t12 - t23) * t55 + (t11 - t166) * t53, t28 * t55 + t116 - t182 + 0.2e1 * t90 + 0.2e1 * t91, -t2 * pkin(4) + t1 * qJ(5) - t11 * t23 + t166 * t12 - t17 * t28 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95 + t180, t5, -t96 ^ 2 - t190, -t12 * t96 - t112;];
tau_reg = t9;
