% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:39
% EndTime: 2019-12-05 17:51:43
% DurationCPUTime: 1.36s
% Computational Cost: add. (2398->237), mult. (4079->317), div. (0->0), fcn. (2468->14), ass. (0->159)
t97 = qJ(1) + pkin(8);
t91 = qJ(3) + t97;
t84 = sin(t91);
t85 = cos(t91);
t197 = g(2) * t85 + g(3) * t84;
t102 = sin(pkin(8));
t192 = pkin(1) * t102;
t151 = qJD(3) * t192;
t104 = cos(pkin(8));
t86 = pkin(1) * t104 + pkin(2);
t198 = qJD(1) * t151 - t86 * qJDD(1);
t106 = sin(qJ(3));
t109 = cos(qJ(3));
t152 = qJD(1) * t192;
t68 = t86 * qJD(1);
t41 = -t106 * t152 + t109 * t68;
t131 = qJD(4) - t41;
t96 = qJD(1) + qJD(3);
t196 = t131 * t96;
t168 = pkin(1) * qJDD(1);
t93 = qJDD(1) + qJDD(3);
t195 = pkin(3) * t93;
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t101 = sin(pkin(9));
t103 = cos(pkin(9));
t59 = -pkin(4) * t103 - pkin(7) * t101 - pkin(3);
t25 = t59 * t96 + t131;
t170 = t106 * t68;
t42 = t109 * t152 + t170;
t33 = qJ(4) * t96 + t42;
t29 = qJD(2) * t101 + t103 * t33;
t7 = -t105 * t29 + t108 * t25;
t171 = t103 * t96;
t71 = -qJD(5) + t171;
t194 = t7 * t71;
t8 = t105 * t25 + t108 * t29;
t193 = t71 * t8;
t191 = t42 * t96;
t161 = qJD(3) * t109;
t46 = -t106 * t151 + t86 * t161;
t45 = qJD(4) + t46;
t190 = t45 * t96;
t51 = t106 * t86 + t109 * t192;
t47 = t51 * qJD(3);
t189 = t47 * t96;
t92 = t96 ^ 2;
t94 = t101 ^ 2;
t188 = t94 * t92;
t160 = qJD(4) * t103;
t163 = t103 * t108;
t164 = t103 * t105;
t43 = -qJ(4) * t164 + t108 * t59;
t166 = qJD(5) * t43;
t187 = -t105 * t42 + t108 * t160 - t41 * t163 + t166;
t44 = qJ(4) * t163 + t105 * t59;
t165 = qJD(5) * t44;
t186 = -t105 * t160 - t108 * t42 + t41 * t164 - t165;
t185 = t197 * t101;
t173 = t103 * t85;
t174 = t103 * t84;
t184 = g(2) * t173 + g(3) * t174;
t182 = g(2) * t84 - g(3) * t85;
t95 = t103 ^ 2;
t181 = t94 + t95;
t98 = t105 ^ 2;
t99 = t108 ^ 2;
t180 = t98 - t99;
t179 = qJ(4) * t93;
t178 = qJD(5) * t8;
t150 = t102 * t168;
t139 = -t198 * t106 + t109 * t150 + t68 * t161;
t18 = qJD(4) * t96 + t139 + t179;
t14 = -t103 * qJDD(2) + t101 * t18;
t177 = t101 * t14;
t28 = -t103 * qJD(2) + t101 * t33;
t176 = t101 * t28;
t175 = t101 * t93;
t172 = t103 * t93;
t67 = -qJDD(5) + t172;
t169 = t67 * t103;
t50 = -t106 * t192 + t109 * t86;
t34 = -t50 + t59;
t48 = qJ(4) + t51;
t20 = t108 * t34 - t48 * t164;
t167 = qJD(5) * t20;
t162 = t105 * t108;
t159 = qJD(5) * t105;
t158 = qJD(5) * t108;
t144 = t101 * t159;
t157 = t7 * t144 + t185;
t156 = pkin(4) * t174;
t155 = pkin(4) * t173;
t138 = qJD(3) * t170 + t106 * t150 + t198 * t109;
t125 = qJDD(4) + t138;
t24 = t125 - t195;
t154 = t24 * t101 - t185;
t153 = t96 * t176;
t149 = t71 * t159;
t148 = t96 * t159;
t147 = t96 * t158;
t146 = -pkin(3) * t84 + t85 * qJ(4);
t145 = t181 * t93;
t143 = -t14 * t103 - g(1);
t142 = t67 - t172;
t141 = t67 + t172;
t140 = (-qJD(5) - t71) * t96;
t137 = t162 * t188;
t135 = t105 * t147;
t107 = sin(qJ(1));
t89 = sin(t97);
t134 = -t107 * pkin(1) - pkin(2) * t89;
t110 = cos(qJ(1));
t90 = cos(t97);
t133 = -t110 * pkin(1) - pkin(2) * t90;
t132 = -t191 - t195;
t130 = -pkin(3) * t85 - t84 * qJ(4);
t129 = g(2) * t110 + g(3) * t107;
t128 = t14 * t48 + t28 * t45;
t127 = t48 * t93 + t190;
t49 = -pkin(3) - t50;
t126 = t49 * t93 + t189;
t15 = qJDD(2) * t101 + t103 * t18;
t124 = t15 * t103 + t177 + t182;
t123 = -t138 + t197;
t122 = -t139 - t182;
t121 = t197 * pkin(7);
t120 = g(1) * t101 - qJD(5) * t25 - t15;
t119 = -qJD(5) * t29 - t153;
t21 = t105 * t34 + t48 * t163;
t118 = t179 + t196;
t117 = t14 * qJ(4) + t131 * t28;
t116 = t134 + t146;
t115 = -t71 ^ 2 - t188;
t13 = t59 * t93 + t125;
t11 = t108 * t13;
t3 = -t105 * t15 + t11 - t178;
t38 = -t85 * t105 + t84 * t163;
t40 = -t105 * t84 - t85 * t163;
t114 = -g(2) * t40 + g(3) * t38 - t103 * t3 + t105 * t177 + t158 * t176;
t113 = t130 + t133;
t2 = qJD(5) * t7 + t105 * t13 + t108 * t15;
t37 = t108 * t85 + t84 * t164;
t39 = -t108 * t84 + t85 * t164;
t112 = -g(2) * t39 - g(3) * t37 + t2 * t103 + t108 * t177 - t144 * t28;
t100 = qJDD(2) - g(1);
t81 = t95 * t93;
t80 = t94 * t93;
t57 = 0.2e1 * t101 * t172;
t54 = t144 * t171;
t36 = (t93 * t99 - 0.2e1 * t135) * t94;
t35 = (t93 * t98 + 0.2e1 * t135) * t94;
t32 = -pkin(3) * t96 + t131;
t27 = 0.2e1 * (t180 * t96 * qJD(5) - t93 * t162) * t94;
t17 = (t141 * t105 + (t71 + t171) * t158) * t101;
t16 = t54 + (-t108 * t141 + t149) * t101;
t5 = -qJD(5) * t21 + t108 * t47 - t45 * t164;
t4 = t105 * t47 + t45 * t163 + t167;
t1 = [0, 0, 0, 0, 0, qJDD(1), t129, -g(2) * t107 + g(3) * t110, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(2) * t90 + g(3) * t89 + 0.2e1 * t104 * t168, -g(2) * t89 + g(3) * t90 - 0.2e1 * t150, 0, (t129 + (t102 ^ 2 + t104 ^ 2) * t168) * pkin(1), 0, 0, 0, 0, 0, t93, t50 * t93 + t123 - t189, -t46 * t96 - t51 * t93 + t122, 0, -g(2) * t133 - g(3) * t134 - t138 * t50 + t139 * t51 - t41 * t47 + t42 * t46, t80, t57, 0, t81, 0, 0, (-t126 - t24) * t103 + t184, t101 * t126 + t154, t48 * t145 + t181 * t190 + t124, t24 * t49 + t32 * t47 - g(2) * t113 - g(3) * t116 + (t15 * t48 + t29 * t45) * t103 + t128 * t101, t36, t27, t16, t35, t17, t169, -t20 * t67 - t5 * t71 + (t105 * t127 + t147 * t48) * t94 + t114, t21 * t67 + t4 * t71 + (t108 * t127 - t148 * t48) * t94 + t112, ((-t21 * t93 - t2 + (-t4 + t167) * t96) * t105 + (-t20 * t93 - t5 * t96 - t3 + (-t21 * t96 - t8) * qJD(5)) * t108) * t101 + t157, t2 * t21 + t8 * t4 + t3 * t20 + t7 * t5 - g(2) * (t113 - t155) - g(3) * (t116 - t156) + (t121 + t128) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101 * t15 + t143, 0, 0, 0, 0, 0, 0, (t142 * t105 + (t71 - t171) * t158) * t101, t54 + (t108 * t142 - t149) * t101, 0, (-t105 * t3 + t108 * t2 + (-t105 * t8 - t108 * t7) * qJD(5)) * t101 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t123 + t191, t41 * t96 + t122, 0, 0, t80, t57, 0, t81, 0, 0, (-t132 - t24) * t103 + t184, t101 * t132 + t154, qJ(4) * t145 + t181 * t196 + t124, -t24 * pkin(3) - t32 * t42 - g(2) * t130 - g(3) * t146 + (t15 * qJ(4) + t131 * t29) * t103 + t117 * t101, t36, t27, t16, t35, t17, t169, -t43 * t67 - t186 * t71 + (qJ(4) * t147 + t105 * t118) * t94 + t114, t44 * t67 + t187 * t71 + (-qJ(4) * t148 + t108 * t118) * t94 + t112, ((-t44 * t93 - t2 + (t166 - t187) * t96) * t105 + (-t178 - t43 * t93 - t3 + (-t165 - t186) * t96) * t108) * t101 + t157, t2 * t44 + t3 * t43 - g(2) * (t130 - t155) - g(3) * (t146 - t156) + t187 * t8 + t186 * t7 + (t121 + t117) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, t175, -t181 * t92, (-t103 * t29 - t176) * t96 + t24 - t197, 0, 0, 0, 0, 0, 0, t105 * t115 - t108 * t67, t105 * t67 + t108 * t115, (-t98 - t99) * t175, -t153 + (t3 - t193) * t108 + (t2 + t194) * t105 - t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, -t180 * t188, (t105 * t140 + t108 * t93) * t101, -t137, (-t105 * t93 + t108 * t140) * t101, -t67, -g(2) * t37 + g(3) * t39 + t105 * t120 + t108 * t119 + t11 - t193, -g(2) * t38 - g(3) * t40 - t194 + t120 * t108 + (-t119 - t13) * t105, 0, 0;];
tau_reg = t1;
