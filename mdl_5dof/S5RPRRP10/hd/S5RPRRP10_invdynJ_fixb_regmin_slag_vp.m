% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:07
% EndTime: 2019-12-31 18:52:12
% DurationCPUTime: 1.91s
% Computational Cost: add. (2456->283), mult. (5853->363), div. (0->0), fcn. (4337->10), ass. (0->150)
t113 = sin(qJ(3));
t116 = cos(qJ(3));
t108 = sin(pkin(8));
t182 = pkin(6) + qJ(2);
t85 = t182 * t108;
t81 = qJD(1) * t85;
t109 = cos(pkin(8));
t86 = t182 * t109;
t82 = qJD(1) * t86;
t50 = -t113 * t81 + t116 * t82;
t208 = t50 * qJD(3);
t107 = pkin(8) + qJ(3);
t100 = cos(t107);
t114 = sin(qJ(1));
t117 = cos(qJ(1));
t141 = g(1) * t117 + g(2) * t114;
t99 = sin(t107);
t122 = -g(3) * t100 + t141 * t99;
t160 = qJD(1) * qJD(2);
t197 = qJDD(1) * t182 + t160;
t61 = t197 * t108;
t62 = t197 * t109;
t131 = -t113 * t62 - t116 * t61 - t208;
t19 = -qJDD(3) * pkin(3) - t131;
t170 = t113 * t108;
t152 = qJD(1) * t170;
t167 = t116 * t109;
t93 = qJD(1) * t167;
t75 = t93 - t152;
t65 = qJD(4) - t75;
t207 = -qJD(4) * pkin(7) * t65 + t122 - t19;
t171 = qJDD(1) * pkin(1);
t202 = g(1) * t114 - g(2) * t117;
t129 = t202 - qJDD(2) + t171;
t206 = t202 * t99;
t112 = sin(qJ(4));
t146 = t112 * t65;
t115 = cos(qJ(4));
t80 = t116 * t108 + t113 * t109;
t76 = t80 * qJD(1);
t58 = t112 * qJD(3) + t115 * t76;
t205 = t58 * t146;
t203 = -t113 * t82 - t116 * t81;
t165 = t117 * t115;
t169 = t114 * t112;
t66 = t100 * t169 + t165;
t166 = t117 * t112;
t168 = t114 * t115;
t68 = -t100 * t166 + t168;
t201 = -g(1) * t68 + g(2) * t66;
t200 = qJ(2) * qJDD(1);
t95 = t109 * pkin(2) + pkin(1);
t84 = -qJD(1) * t95 + qJD(2);
t28 = -t75 * pkin(3) - t76 * pkin(7) + t84;
t44 = qJD(3) * pkin(7) + t50;
t13 = t112 * t28 + t115 * t44;
t158 = t109 * qJDD(1);
t159 = t108 * qJDD(1);
t156 = qJD(3) * t93 + t113 * t158 + t116 * t159;
t125 = qJD(3) * t152 - t156;
t137 = t113 * t159 - t116 * t158;
t78 = t80 * qJD(3);
t48 = qJD(1) * t78 + t137;
t83 = -qJDD(1) * t95 + qJDD(2);
t20 = t48 * pkin(3) + pkin(7) * t125 + t83;
t17 = t115 * t20;
t134 = -t113 * t61 + t116 * t62;
t18 = qJDD(3) * pkin(7) + qJD(3) * t203 + t134;
t161 = t115 * qJD(3);
t163 = qJD(4) * t112;
t22 = -qJD(4) * t161 - t112 * qJDD(3) + t115 * t125 + t163 * t76;
t42 = qJDD(4) + t48;
t1 = t42 * pkin(4) + t22 * qJ(5) - qJD(4) * t13 - t58 * qJD(5) - t112 * t18 + t17;
t56 = t112 * t76 - t161;
t8 = -t56 * qJ(5) + t13;
t199 = t65 * t8 + t1;
t192 = g(3) * t99;
t198 = t141 * t100 + t192;
t196 = t58 ^ 2;
t12 = -t112 * t44 + t115 * t28;
t7 = -t58 * qJ(5) + t12;
t6 = t65 * pkin(4) + t7;
t195 = -t7 + t6;
t181 = qJ(5) + pkin(7);
t149 = qJD(4) * t181;
t172 = qJ(5) * t115;
t45 = t76 * pkin(3) - t75 * pkin(7);
t35 = t115 * t45;
t191 = -t76 * pkin(4) - t115 * t149 + t172 * t75 - t35 + (-qJD(5) + t203) * t112;
t190 = pkin(4) * t112;
t185 = t56 * t75;
t184 = t58 * t76;
t183 = t76 * t56;
t162 = qJD(4) * t115;
t121 = -t115 * qJDD(3) - t112 * t125;
t23 = qJD(4) * t58 + t121;
t180 = -t112 * t23 - t56 * t162;
t179 = t112 * t45 + t115 * t203;
t79 = -t167 + t170;
t47 = t79 * pkin(3) - t80 * pkin(7) - t95;
t54 = -t113 * t85 + t116 * t86;
t51 = t115 * t54;
t178 = t112 * t47 + t51;
t173 = qJ(5) * t112;
t177 = t115 * qJD(5) - t112 * t149 + t173 * t75 - t179;
t176 = t112 * t42;
t174 = t22 * t112;
t164 = t108 ^ 2 + t109 ^ 2;
t53 = t113 * t86 + t116 * t85;
t29 = -qJD(2) * t79 - qJD(3) * t53;
t77 = t79 * qJD(3);
t46 = t78 * pkin(3) + t77 * pkin(7);
t157 = t112 * t46 + t115 * t29 + t47 * t162;
t154 = t80 * t163;
t153 = t80 * t162;
t151 = t182 + t190;
t148 = -qJD(4) * t28 - t18;
t147 = t115 * t65;
t145 = t164 * qJD(1) ^ 2;
t144 = 0.2e1 * t164;
t143 = -t162 * t44 + t17;
t126 = t112 * t20 + t115 * t18 + t28 * t162 - t163 * t44;
t2 = -t23 * qJ(5) - t56 * qJD(5) + t126;
t142 = -t6 * t65 + t2;
t43 = -qJD(3) * pkin(3) - t203;
t139 = t19 * t80 - t43 * t77;
t138 = t42 * t80 - t65 * t77;
t97 = t115 * pkin(4) + pkin(3);
t136 = t100 * t97 + t181 * t99;
t132 = qJ(5) * t77 - qJD(5) * t80;
t130 = t115 * t42 + t75 * t146 - t163 * t65;
t128 = t136 + t95;
t127 = -pkin(7) * t42 + t43 * t65;
t124 = t129 + t171;
t5 = t23 * pkin(4) + qJDD(5) + t19;
t120 = t144 * t160 - t141;
t30 = qJD(2) * t80 + qJD(3) * t54;
t88 = t181 * t115;
t87 = t181 * t112;
t69 = t100 * t165 + t169;
t67 = -t100 * t168 + t166;
t55 = t56 ^ 2;
t38 = t115 * t47;
t36 = t115 * t46;
t24 = t56 * pkin(4) + qJD(5) + t43;
t14 = -t173 * t80 + t178;
t10 = t79 * pkin(4) - t112 * t54 - t172 * t80 + t38;
t4 = -qJ(5) * t153 + (-qJD(4) * t54 + t132) * t112 + t157;
t3 = t78 * pkin(4) - t112 * t29 + t36 + t132 * t115 + (-t51 + (qJ(5) * t80 - t47) * t112) * qJD(4);
t9 = [qJDD(1), t202, t141, t124 * t109, -t124 * t108, t144 * t200 + t120, pkin(1) * t129 + (t164 * t200 + t120) * qJ(2), -t125 * t80 - t76 * t77, t125 * t79 - t80 * t48 - t77 * t75 - t76 * t78, -t77 * qJD(3) + t80 * qJDD(3), -t78 * qJD(3) - t79 * qJDD(3), 0, -t30 * qJD(3) - t53 * qJDD(3) + t100 * t202 - t95 * t48 + t84 * t78 + t83 * t79, -t29 * qJD(3) - t54 * qJDD(3) + t125 * t95 - t84 * t77 + t83 * t80 - t206, -t58 * t154 + (-t22 * t80 - t58 * t77) * t115, -(-t112 * t58 - t115 * t56) * t77 + (t174 - t115 * t23 + (t112 * t56 - t115 * t58) * qJD(4)) * t80, t115 * t138 - t154 * t65 - t22 * t79 + t58 * t78, -t112 * t138 - t153 * t65 - t23 * t79 - t56 * t78, t42 * t79 + t65 * t78, (-t162 * t54 + t36) * t65 + t38 * t42 + t143 * t79 + t12 * t78 + t30 * t56 + t53 * t23 + t43 * t153 - g(1) * t67 - g(2) * t69 + ((-qJD(4) * t47 - t29) * t65 - t54 * t42 + t148 * t79 + t139) * t112, -(-t163 * t54 + t157) * t65 - t178 * t42 - t126 * t79 - t13 * t78 + t30 * t58 - t53 * t22 - t43 * t154 - g(1) * t66 - g(2) * t68 + t139 * t115, t10 * t22 - t14 * t23 - t3 * t58 - t4 * t56 + t206 - (-t112 * t8 - t115 * t6) * t77 + (-t1 * t115 - t2 * t112 + (t112 * t6 - t115 * t8) * qJD(4)) * t80, t2 * t14 + t8 * t4 + t1 * t10 + t6 * t3 + t5 * (t190 * t80 + t53) + t24 * ((-t112 * t77 + t153) * pkin(4) + t30) + (-g(1) * t151 - g(2) * t128) * t117 + (g(1) * t128 - g(2) * t151) * t114; 0, 0, 0, -t158, t159, -t145, -qJ(2) * t145 - t129, 0, 0, 0, 0, 0, 0.2e1 * t76 * qJD(3) + t137, (t75 - t152) * qJD(3) + t156, 0, 0, 0, 0, 0, t130 - t183, -t115 * t65 ^ 2 - t176 - t184, (t22 + t185) * t115 + t205 + t180, t142 * t112 + t199 * t115 - t24 * t76 - t202; 0, 0, 0, 0, 0, 0, 0, -t76 * t75, -t75 ^ 2 + t76 ^ 2, (-t75 - t152) * qJD(3) + t156, -t137, qJDD(3), -t84 * t76 + t122 + t131 + t208, -t84 * t75 - t134 + t198, t147 * t58 - t174, (-t22 + t185) * t115 - t205 + t180, t147 * t65 + t176 - t184, t130 + t183, -t65 * t76, -pkin(3) * t23 - t12 * t76 - t35 * t65 - t50 * t56 + (t203 * t65 + t127) * t112 + t207 * t115, pkin(3) * t22 - t207 * t112 + t127 * t115 + t13 * t76 + t179 * t65 - t50 * t58, -t199 * t112 + t142 * t115 - t177 * t56 - t191 * t58 - t87 * t22 - t88 * t23 - t198, t2 * t88 - t1 * t87 - t5 * t97 - g(3) * t136 + t177 * t8 + t191 * t6 + (pkin(4) * t146 - t50) * t24 + t141 * (-t100 * t181 + t97 * t99); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t56, -t55 + t196, t56 * t65 - t22, -t121 + (-qJD(4) + t65) * t58, t42, t13 * t65 - t43 * t58 + (t148 + t192) * t112 + t143 + t201, g(1) * t69 - g(2) * t67 + t115 * t192 + t12 * t65 + t43 * t56 - t126, pkin(4) * t22 - t195 * t56, t195 * t8 + (t112 * t192 - t24 * t58 + t1 + t201) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 - t196, t8 * t56 + t6 * t58 - t122 + t5;];
tau_reg = t9;
