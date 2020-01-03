% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:34
% EndTime: 2019-12-31 18:47:37
% DurationCPUTime: 1.27s
% Computational Cost: add. (1881->236), mult. (2788->245), div. (0->0), fcn. (1437->6), ass. (0->140)
t87 = qJDD(1) - qJDD(3);
t93 = cos(qJ(4));
t90 = t93 ^ 2;
t168 = t90 * t87;
t91 = sin(qJ(4));
t89 = t91 ^ 2;
t170 = t87 * t89;
t206 = -t168 - t170;
t179 = sin(qJ(1));
t180 = cos(qJ(1));
t92 = sin(qJ(3));
t94 = cos(qJ(3));
t47 = -t179 * t92 - t180 * t94;
t48 = -t179 * t94 + t180 * t92;
t123 = g(1) * t47 + g(2) * t48;
t141 = qJD(1) - qJD(3);
t205 = t141 ^ 2;
t158 = t89 + t90;
t121 = pkin(4) * t93 + qJ(5) * t91;
t117 = pkin(3) + t121;
t204 = t117 * t87;
t199 = t117 * t141;
t148 = qJ(2) * qJD(1);
t95 = -pkin(1) - pkin(2);
t61 = t95 * qJD(1) + qJD(2);
t34 = -t92 * t148 + t94 * t61;
t13 = -t34 + t199;
t203 = t13 * t141;
t188 = pkin(3) * t141;
t24 = -t34 + t188;
t202 = t141 * t24;
t201 = t141 * t34;
t35 = t94 * t148 + t61 * t92;
t172 = t35 * t141;
t200 = -qJD(3) * t148 + t95 * qJDD(1) + qJDD(2);
t198 = qJD(4) * t141;
t143 = qJD(1) * qJD(2);
t144 = qJ(2) * qJDD(1);
t197 = qJD(3) * t61 + t143 + t144;
t11 = t197 * t94 + t200 * t92;
t9 = -pkin(7) * t87 + t11;
t138 = t158 * t9;
t55 = t94 * qJ(2) + t92 * t95;
t196 = t121 * t47;
t195 = t121 * t48;
t54 = -t92 * qJ(2) + t94 * t95;
t131 = -qJD(4) * pkin(4) + qJD(5);
t25 = -pkin(7) * t141 + t35;
t176 = t25 * t91;
t15 = t131 + t176;
t147 = qJD(4) * qJ(5);
t175 = t25 * t93;
t16 = t147 + t175;
t118 = t15 * t91 + t16 * t93;
t142 = qJDD(4) * qJ(5);
t8 = t93 * t9;
t2 = t142 + t8 + (qJD(5) - t176) * qJD(4);
t150 = qJDD(4) * pkin(4);
t152 = qJD(4) * t93;
t7 = t91 * t9;
t3 = t25 * t152 + qJDD(5) - t150 + t7;
t122 = t2 * t93 + t3 * t91;
t194 = g(3) * t93 - t123 * t91;
t136 = t158 * t34;
t193 = t206 * pkin(7) + t136 * t141;
t192 = -qJDD(4) * t92 + 0.2e1 * t94 * t198;
t102 = (t15 * t93 - t16 * t91) * qJD(4) + t122;
t120 = pkin(4) * t91 - qJ(5) * t93;
t36 = t120 * qJD(4) - t91 * qJD(5);
t50 = -pkin(7) + t55;
t146 = qJDD(4) * t50;
t27 = t117 - t54;
t32 = t94 * qJD(2) + t54 * qJD(3);
t191 = (-t141 * t27 - t13 - t32) * qJD(4) - t146;
t159 = -t89 + t90;
t166 = t93 * t87;
t190 = -0.2e1 * t159 * t198 - 0.2e1 * t91 * t166;
t189 = pkin(3) * t87;
t182 = t47 * pkin(3);
t181 = t48 * pkin(3);
t174 = t32 * t141;
t33 = t92 * qJD(2) + t55 * qJD(3);
t173 = t33 * t141;
t171 = t36 * t141;
t169 = t141 * t91;
t167 = t92 * t87;
t165 = t94 * t87;
t153 = qJD(4) * t91;
t163 = t34 * t153 - t93 * t172;
t162 = t36 - t35;
t161 = t180 * pkin(1) + t179 * qJ(2);
t160 = g(1) * t179 - g(2) * t180;
t156 = pkin(1) * qJDD(1);
t155 = pkin(7) * qJDD(4);
t140 = t180 * pkin(2) + t161;
t139 = 0.2e1 * t143;
t137 = t24 + t188;
t134 = t13 + t199;
t132 = -t7 + t194;
t129 = qJDD(2) - t156;
t126 = -0.2e1 * t152 * t169;
t125 = -t179 * pkin(1) + t180 * qJ(2);
t12 = -t197 * t92 + t200 * t94;
t124 = g(1) * t48 - g(2) * t47;
t116 = t48 * pkin(7) + t140 - t182;
t96 = qJD(4) ^ 2;
t115 = pkin(7) * t96 - t124;
t114 = -t50 * t96 - t124;
t113 = g(1) * t180 + g(2) * t179;
t112 = -t179 * pkin(2) + t125;
t111 = -t158 * t174 + t206 * t50 - t123;
t110 = t138 + t123;
t109 = t141 * t158;
t108 = t124 + t12;
t49 = pkin(3) - t54;
t107 = -t146 + (-t141 * t49 - t24 - t32) * qJD(4);
t106 = -t165 + (-t96 - t205) * t92;
t10 = -t12 + t189;
t105 = -t10 - t115 - t189;
t1 = -t12 - t171 + t204;
t104 = -t1 - t115 - t204;
t103 = t47 * pkin(7) + t112 + t181;
t14 = t33 - t36;
t101 = t14 * t141 + t27 * t87 + t1 + t114;
t100 = -t49 * t87 - t10 - t114 - t173;
t99 = -t11 - t123;
t98 = t102 + t123;
t97 = qJD(1) ^ 2;
t65 = t91 * t87;
t59 = t91 * t205 * t93;
t58 = qJDD(4) * t93 - t96 * t91;
t57 = qJDD(4) * t91 + t93 * t96;
t40 = t159 * t205;
t39 = t120 * t141;
t31 = t126 + t168;
t30 = t126 - t170;
t6 = t109 * t141 * t94 - t158 * t167;
t5 = t106 * t93 + t192 * t91;
t4 = t106 * t91 - t192 * t93;
t17 = [0, 0, 0, 0, 0, qJDD(1), t160, t113, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t156 + t160, 0, -t113 + t139 + 0.2e1 * t144, -t129 * pkin(1) - g(1) * t125 - g(2) * t161 + (t139 + t144) * qJ(2), 0, 0, 0, 0, 0, t87, -t54 * t87 - t108 + t173, t55 * t87 + t174 - t99, 0, -g(1) * t112 - g(2) * t140 + t11 * t55 + t12 * t54 + t35 * t32 - t34 * t33, -t30, -t190, -t57, t31, -t58, 0, -t100 * t93 + t107 * t91, t100 * t91 + t107 * t93, t111 - t138, t158 * t32 * t25 - g(1) * t103 - g(2) * t116 + t10 * t49 + t50 * t138 + t24 * t33, -t30, -t57, t190, 0, t58, t31, t101 * t93 + t191 * t91, -t102 + t111, t101 * t91 - t191 * t93, t1 * t27 + t13 * t14 - g(1) * (t103 + t195) - g(2) * (t116 - t196) + t118 * t32 + (t15 * t152 - t16 * t153 + t122) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t97, -qJ(2) * t97 + t129 - t160, 0, 0, 0, 0, 0, 0, -t205 * t92 - t165, -t205 * t94 + t167, 0, (t12 - t172) * t94 + (t11 + t201) * t92 - t160, 0, 0, 0, 0, 0, 0, t5, -t4, t6, (t138 - t202) * t92 + (-t109 * t25 - t10) * t94 - t160, 0, 0, 0, 0, 0, 0, t5, t6, t4, (-t141 * t118 - t1) * t94 + (t102 - t203) * t92 - t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t108 - t172, t99 - t201, 0, 0, t30, t190, t57, -t31, t58, 0, (t137 * qJD(4) - t155) * t91 + t105 * t93 + t163, (-t155 + (t137 + t34) * qJD(4)) * t93 + (-t105 + t172) * t91, t110 + t193, -t10 * pkin(3) + t110 * pkin(7) + g(1) * t181 - g(2) * t182 - t25 * t136 - t24 * t35, t30, t57, -t190, 0, -t58, -t31, (t134 * qJD(4) - t155) * t91 + (t104 + t171) * t93 + t163, t98 + t193, (t155 + (-t134 - t34) * qJD(4)) * t93 + (t141 * t162 + t104) * t91, -t1 * t117 - g(1) * (-t181 - t195) - g(2) * (t182 + t196) - t118 * t34 + t162 * t13 + t98 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t40, -t65, t59, -t166, qJDD(4), t24 * t169 + t132, -g(3) * t91 - t8 + (-t123 + t202) * t93, 0, 0, -t59, -t65, t40, qJDD(4), t166, t59, 0.2e1 * t150 - qJDD(5) - (-t13 * t91 - t39 * t93) * t141 + t132, t120 * t87 - ((t16 - t147) * t91 + (t131 - t15) * t93) * t141, 0.2e1 * t142 + 0.2e1 * qJD(4) * qJD(5) + t8 + (t141 * t39 + g(3)) * t91 + (t123 - t203) * t93, t2 * qJ(5) - t3 * pkin(4) + t13 * t39 - t15 * t175 + g(3) * t121 + (qJD(5) + t176) * t16 - t123 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t59, -t65, -t205 * t89 - t96, -qJD(4) * t16 - t13 * t169 - t194 + t3;];
tau_reg = t17;
