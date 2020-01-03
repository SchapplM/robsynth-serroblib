% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:08
% EndTime: 2019-12-31 17:17:11
% DurationCPUTime: 1.17s
% Computational Cost: add. (1593->243), mult. (3752->306), div. (0->0), fcn. (2437->8), ass. (0->131)
t92 = qJD(2) + qJD(3);
t171 = pkin(6) + pkin(5);
t95 = qJ(2) + qJ(3);
t89 = sin(t95);
t90 = cos(t95);
t176 = t90 * pkin(3) + t89 * qJ(4);
t91 = qJDD(2) + qJDD(3);
t86 = t91 * qJ(4);
t87 = t92 * qJD(4);
t175 = t86 + t87;
t96 = sin(qJ(3));
t99 = cos(qJ(2));
t159 = t96 * t99;
t168 = cos(qJ(3));
t97 = sin(qJ(2));
t53 = t168 * t97 + t159;
t47 = t53 * qJD(1);
t88 = t91 * pkin(3);
t174 = qJDD(4) - t88;
t100 = cos(qJ(1));
t98 = sin(qJ(1));
t127 = g(1) * t98 - g(2) * t100;
t59 = t171 * t97;
t60 = t171 * t99;
t119 = -t168 * t59 - t96 * t60;
t139 = qJD(2) * t171;
t56 = t97 * t139;
t13 = t119 * qJD(3) - t139 * t159 - t168 * t56;
t39 = t168 * t60 - t96 * t59;
t173 = t127 * t89 + t13 * t92 + t39 * t91;
t172 = t47 ^ 2;
t170 = g(3) * t99;
t169 = t99 * pkin(2);
t57 = qJD(1) * t60;
t161 = t96 * t57;
t154 = qJD(2) * pkin(2);
t55 = qJD(1) * t59;
t51 = -t55 + t154;
t26 = t168 * t51 - t161;
t167 = t26 * t92;
t141 = t168 * t57;
t27 = t96 * t51 + t141;
t166 = t27 * t92;
t140 = t168 * t99;
t130 = qJD(1) * t140;
t150 = qJD(1) * t97;
t142 = t96 * t150;
t45 = -t130 + t142;
t165 = t47 * t45;
t163 = t89 * t98;
t162 = t90 * t98;
t160 = t96 * t97;
t136 = t168 * qJD(3);
t31 = -t168 * t55 - t161;
t158 = pkin(2) * t136 + qJD(4) - t31;
t93 = t97 ^ 2;
t94 = t99 ^ 2;
t156 = t93 - t94;
t155 = t93 + t94;
t153 = t100 * t89;
t152 = t100 * t90;
t151 = pkin(5) * qJDD(1);
t149 = qJD(3) * t96;
t148 = qJD(4) - t26;
t147 = t97 * qJDD(1);
t146 = t99 * qJDD(1);
t145 = qJD(1) * qJD(2);
t144 = t97 * t154;
t103 = qJD(1) ^ 2;
t143 = t97 * t103 * t99;
t18 = t45 ^ 2 - t172;
t85 = pkin(1) + t169;
t138 = t97 * t145;
t137 = t99 * t145;
t135 = qJDD(1) * t168;
t36 = qJDD(2) * pkin(2) + t171 * (-t137 - t147);
t37 = t171 * (-t138 + t146);
t5 = t51 * t136 - t57 * t149 + t168 * t37 + t96 * t36;
t134 = t57 * t136 + t51 * t149 - t168 * t36 + t96 * t37;
t133 = -t92 * t130 - t97 * t135 - t96 * t146;
t132 = t97 * t137;
t131 = qJD(2) * t140;
t30 = -t96 * t55 + t141;
t129 = pkin(2) * t149 - t30;
t128 = -pkin(2) * t97 - pkin(3) * t89;
t126 = g(1) * t100 + g(2) * t98;
t125 = -t99 * t135 + t96 * t147;
t58 = t85 * qJD(1);
t124 = t92 * t160;
t24 = t47 * pkin(3) + t45 * qJ(4);
t33 = t92 * t53;
t17 = t33 * qJD(1) + t125;
t52 = -t140 + t160;
t121 = t17 * t52 + t45 * t33;
t120 = t92 * t33 + t91 * t52;
t118 = g(1) * t152 + g(2) * t162 + g(3) * t89 - t5;
t42 = pkin(2) * t138 - t85 * qJDD(1);
t117 = -0.2e1 * pkin(1) * t145 - pkin(5) * qJDD(2);
t114 = g(1) * t153 + g(2) * t163 - g(3) * t90 - t134;
t20 = t45 * pkin(3) - t47 * qJ(4) - t58;
t113 = -t20 * t45 - t118;
t112 = -t58 * t45 + t118;
t14 = t39 * qJD(3) + t171 * t131 - t96 * t56;
t111 = g(1) * t162 - g(2) * t152 + t119 * t91 - t14 * t92;
t16 = qJD(1) * t124 + t133;
t32 = -t99 * t136 + t124 - t131;
t110 = t16 * t52 - t53 * t17 + t32 * t45 - t47 * t33;
t102 = qJD(2) ^ 2;
t109 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t102 + t127;
t108 = pkin(1) * t103 + t126 - t151;
t107 = t58 * t47 + t114;
t106 = t119 * t16 - t13 * t45 + t14 * t47 - t39 * t17 - t126;
t105 = t20 * t47 - t114 + t174;
t84 = -t168 * pkin(2) - pkin(3);
t80 = t96 * pkin(2) + qJ(4);
t63 = t100 * t85;
t62 = qJ(4) * t152;
t61 = qJ(4) * t162;
t25 = t52 * pkin(3) - t53 * qJ(4) - t85;
t23 = t92 * qJ(4) + t27;
t22 = pkin(2) * t150 + t24;
t21 = -t92 * pkin(3) + t148;
t15 = -t32 * t92 + t53 * t91;
t11 = -t47 * t92 + t17;
t10 = -t133 + (-t142 + t45) * t92;
t7 = t33 * pkin(3) + t32 * qJ(4) - t53 * qJD(4) + t144;
t4 = t134 + t174;
t3 = t5 + t175;
t2 = -t16 * t53 - t47 * t32;
t1 = t17 * pkin(3) + t16 * qJ(4) - t47 * qJD(4) + t42;
t6 = [0, 0, 0, 0, 0, qJDD(1), t127, t126, 0, 0, t93 * qJDD(1) + 0.2e1 * t132, -0.2e1 * t156 * t145 + 0.2e1 * t97 * t146, qJDD(2) * t97 + t102 * t99, t94 * qJDD(1) - 0.2e1 * t132, qJDD(2) * t99 - t102 * t97, 0, t109 * t99 + t117 * t97, -t109 * t97 + t117 * t99, 0.2e1 * t151 * t155 - t126, -g(1) * (-t98 * pkin(1) + t100 * pkin(5)) - g(2) * (t100 * pkin(1) + t98 * pkin(5)) + (t155 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), t2, t110, t15, t121, -t120, 0, t45 * t144 - t85 * t17 - t58 * t33 + t42 * t52 + t111, t47 * t144 + t85 * t16 + t58 * t32 + t42 * t53 - t173, t134 * t53 + t26 * t32 - t27 * t33 - t5 * t52 + t106, t5 * t39 + t27 * t13 - t134 * t119 - t26 * t14 - t42 * t85 - t58 * t144 - g(1) * (t100 * t171 - t98 * t85) - g(2) * (t171 * t98 + t63), t2, t15, -t110, 0, t120, t121, t1 * t52 + t25 * t17 + t20 * t33 + t7 * t45 + t111, -t21 * t32 - t23 * t33 - t3 * t52 + t4 * t53 + t106, -t1 * t53 + t25 * t16 + t20 * t32 - t7 * t47 + t173, -g(2) * t63 + t1 * t25 + t23 * t13 + t21 * t14 + t20 * t7 + t3 * t39 - t4 * t119 + (-g(1) * t171 - g(2) * t176) * t100 + (-g(1) * (-t176 - t85) - g(2) * t171) * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t156 * t103, t147, t143, t146, qJDD(2), t108 * t97 - t170, g(3) * t97 + t108 * t99, 0, 0, t165, -t18, t10, -t165, -t125, t91, t30 * t92 + (-t92 * t149 - t45 * t150 + t168 * t91) * pkin(2) + t107, t31 * t92 + (-t92 * t136 - t47 * t150 - t91 * t96) * pkin(2) + t112, (t27 - t30) * t47 + (-t26 + t31) * t45 + (t168 * t16 - t17 * t96 + (-t168 * t45 + t47 * t96) * qJD(3)) * pkin(2), t26 * t30 - t27 * t31 + (-t168 * t134 - t170 + t5 * t96 + (t168 * t27 - t26 * t96) * qJD(3) + (qJD(1) * t58 + t126) * t97) * pkin(2), t165, t10, t18, t91, t11, -t165, -t129 * t92 - t22 * t45 - t84 * t91 - t105, -t84 * t16 - t80 * t17 + (t129 + t23) * t47 + (t21 - t158) * t45, t158 * t92 + t22 * t47 + t80 * t91 + t113 + t175, t3 * t80 + t4 * t84 - t20 * t22 - g(1) * (t100 * t128 + t62) - g(2) * (t128 * t98 + t61) - g(3) * (t176 + t169) + t158 * t23 + t129 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, -t18, t10, -t165, -t125, t91, t107 + t166, t112 + t167, 0, 0, t165, t10, t18, t91, t11, -t165, -t24 * t45 - t105 + t166 + t88, pkin(3) * t16 - t17 * qJ(4) + (t23 - t27) * t47 + (t21 - t148) * t45, t24 * t47 + t113 - t167 + 0.2e1 * t86 + 0.2e1 * t87, t3 * qJ(4) - t4 * pkin(3) - t20 * t24 - t21 * t27 - g(1) * (-pkin(3) * t153 + t62) - g(2) * (-pkin(3) * t163 + t61) - g(3) * t176 + t148 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91 + t165, t10, -t92 ^ 2 - t172, -t23 * t92 + t105;];
tau_reg = t6;
