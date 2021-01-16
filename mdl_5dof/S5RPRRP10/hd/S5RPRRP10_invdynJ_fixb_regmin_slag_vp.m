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
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 19:14:49
% EndTime: 2021-01-15 19:15:06
% DurationCPUTime: 2.78s
% Computational Cost: add. (3258->330), mult. (7740->405), div. (0->0), fcn. (5766->10), ass. (0->171)
t234 = 2 * qJD(3);
t117 = cos(pkin(8));
t124 = cos(qJ(3));
t116 = sin(pkin(8));
t121 = sin(qJ(3));
t182 = t116 * t121;
t230 = t124 * t117 - t182;
t233 = t230 * qJD(1);
t119 = pkin(6) + qJ(2);
t93 = t119 * t116;
t87 = qJD(1) * t93;
t94 = t119 * t117;
t88 = qJD(1) * t94;
t56 = -t121 * t87 + t124 * t88;
t232 = qJD(3) * t56;
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t86 = t116 * t124 + t117 * t121;
t135 = t86 * qJDD(1);
t173 = qJD(4) * t123;
t175 = qJD(3) * t120;
t82 = t86 * qJD(1);
t26 = -t123 * qJDD(3) + t120 * t135 + t82 * t173 + (qJD(4) + t233) * t175;
t64 = t123 * t82 + t175;
t231 = t120 * t64;
t183 = qJDD(1) * pkin(1);
t122 = sin(qJ(1));
t125 = cos(qJ(1));
t223 = g(1) * t122 - g(2) * t125;
t140 = -qJDD(2) + t183 + t223;
t103 = t117 * pkin(2) + pkin(1);
t92 = -t103 * qJD(1) + qJD(2);
t32 = -pkin(3) * t233 - pkin(7) * t82 + t92;
t50 = qJD(3) * pkin(7) + t56;
t15 = t120 * t32 + t123 * t50;
t83 = t230 * qJD(3);
t134 = qJD(1) * t83;
t224 = -pkin(7) * t86 - t103;
t170 = t117 * qJDD(1);
t148 = qJDD(1) * t182 - t124 * t170;
t84 = t86 * qJD(3);
t54 = qJD(1) * t84 + t148;
t23 = t54 * pkin(3) - pkin(7) * t134 + t224 * qJDD(1) + qJDD(2);
t19 = t123 * t23;
t171 = qJD(1) * qJD(2);
t219 = t119 * qJDD(1) + t171;
t67 = t219 * t116;
t68 = t219 * t117;
t147 = -t121 * t67 + t124 * t68;
t226 = -t121 * t88 - t124 * t87;
t20 = qJDD(3) * pkin(7) + qJD(3) * t226 + t147;
t132 = -qJD(4) * t15 - t120 * t20 + t19;
t128 = t135 + t134;
t172 = t123 * qJD(3);
t174 = qJD(4) * t120;
t25 = -qJD(4) * t172 - t120 * qJDD(3) - t123 * t128 + t82 * t174;
t192 = qJ(5) * t25;
t48 = qJDD(4) + t54;
t216 = pkin(4) * t48;
t1 = -qJD(5) * t64 + t132 + t192 + t216;
t62 = t120 * t82 - t172;
t10 = -qJ(5) * t62 + t15;
t71 = qJD(4) - t233;
t204 = t10 * t71;
t228 = t1 + t204;
t227 = t71 * t231;
t115 = pkin(8) + qJ(3);
t107 = sin(t115);
t225 = t223 * t107;
t152 = g(1) * t125 + g(2) * t122;
t139 = t152 * t107;
t222 = qJ(2) * qJDD(1);
t205 = g(3) * t120;
t108 = cos(t115);
t178 = t123 * t125;
t181 = t120 * t122;
t72 = t108 * t181 + t178;
t179 = t122 * t123;
t180 = t120 * t125;
t74 = -t108 * t180 + t179;
t221 = -g(1) * t74 + g(2) * t72 + t107 * t205;
t207 = g(3) * t107;
t220 = t152 * t108 + t207;
t218 = t64 ^ 2;
t14 = -t120 * t50 + t123 * t32;
t9 = -qJ(5) * t64 + t14;
t8 = pkin(4) * t71 + t9;
t217 = -t9 + t8;
t215 = pkin(4) * t62;
t211 = pkin(4) * t120;
t206 = g(3) * t108;
t203 = t62 * t71;
t202 = t62 * t233;
t201 = t62 * t82;
t200 = t64 * t71;
t199 = t64 * t82;
t118 = -qJ(5) - pkin(7);
t198 = -t120 * t26 - t62 * t173;
t51 = pkin(3) * t82 - pkin(7) * t233;
t197 = t120 * t51 + t123 * t226;
t53 = -pkin(3) * t230 + t224;
t60 = -t121 * t93 + t124 * t94;
t57 = t123 * t60;
t196 = t120 * t53 + t57;
t162 = qJD(4) * t118;
t185 = qJ(5) * t120;
t195 = qJD(5) * t123 + t120 * t162 + t185 * t233 - t197;
t184 = qJ(5) * t123;
t41 = t123 * t51;
t194 = -pkin(4) * t82 + t123 * t162 + t233 * t184 - t41 + (-qJD(5) + t226) * t120;
t193 = (g(1) * t178 + g(2) * t179) * t107;
t191 = qJ(5) * t26;
t190 = t120 * t25;
t189 = t120 * t48;
t187 = t120 * t233;
t176 = t116 ^ 2 + t117 ^ 2;
t59 = t121 * t94 + t124 * t93;
t33 = qJD(2) * t230 - qJD(3) * t59;
t52 = pkin(3) * t84 - pkin(7) * t83;
t169 = t120 * t52 + t123 * t33 + t53 * t173;
t168 = t86 * t174;
t167 = t86 * t173;
t144 = -t121 * t68 - t124 * t67 - t232;
t21 = -qJDD(3) * pkin(3) - t144;
t5 = pkin(4) * t26 + qJDD(5) + t21;
t166 = -t5 - t206;
t165 = -qJD(5) - t215;
t164 = t119 + t211;
t161 = t120 * t23 + t123 * t20 + t32 * t173 - t50 * t174;
t160 = t123 * t71;
t159 = t176 * qJD(1) ^ 2;
t158 = 0.2e1 * t176;
t157 = qJD(4) * pkin(7) * t71 + t21;
t156 = -g(1) * t72 - g(2) * t74;
t73 = -t108 * t179 + t180;
t75 = t108 * t178 + t181;
t155 = -g(1) * t73 - g(2) * t75;
t49 = -qJD(3) * pkin(3) - t226;
t27 = -t165 + t49;
t154 = t27 * t83 + t5 * t86;
t2 = -qJD(5) * t62 + t161 - t191;
t153 = -t71 * t8 + t2;
t150 = t21 * t86 + t49 * t83;
t149 = t48 * t86 + t71 * t83;
t145 = -qJ(5) * t83 - qJD(5) * t86;
t105 = pkin(4) * t123 + pkin(3);
t143 = t105 * t108 - t107 * t118;
t141 = t123 * t48 + (-t174 + t187) * t71;
t138 = -pkin(7) * t48 + t71 * t49;
t137 = t103 + t143;
t133 = g(1) * t75 - g(2) * t73 + t123 * t207 - t161;
t131 = t158 * t171 - t152;
t34 = qJD(2) * t86 + qJD(3) * t60;
t129 = t132 + t221;
t99 = t108 * t205;
t96 = t118 * t123;
t95 = t118 * t120;
t91 = -t103 * qJDD(1) + qJDD(2);
t61 = t62 ^ 2;
t44 = t123 * t53;
t42 = t123 * t52;
t35 = t86 * t211 + t59;
t28 = pkin(4) * t187 + t56;
t22 = (t120 * t83 + t167) * pkin(4) + t34;
t16 = -t86 * t185 + t196;
t12 = -pkin(4) * t230 - t120 * t60 - t86 * t184 + t44;
t7 = -t71 ^ 2 * t123 - t189 - t199;
t6 = t141 - t201;
t4 = -qJ(5) * t167 + (-qJD(4) * t60 + t145) * t120 + t169;
t3 = pkin(4) * t84 - t120 * t33 + t42 + t145 * t123 + (-t57 + (qJ(5) * t86 - t53) * t120) * qJD(4);
t11 = [qJDD(1), t223, t152, (t140 + t183) * t117, t158 * t222 + t131, t140 * pkin(1) + (t176 * t222 + t131) * qJ(2), t128 * t86 + t82 * t83, t128 * t230 + t233 * t83 - t86 * t54 - t82 * t84, qJD(3) * t83 + qJDD(3) * t86, -qJD(3) * t84 + qJDD(3) * t230, 0, -qJD(3) * t34 - qJDD(3) * t59 - t103 * t54 + t108 * t223 - t230 * t91 + t84 * t92, -t60 * qJDD(3) + t92 * t83 + t91 * t86 - t225 - t103 * t135 + (-t103 * t233 - t33) * qJD(3), -t64 * t168 + (-t25 * t86 + t64 * t83) * t123, (-t123 * t62 - t231) * t83 + (t190 - t123 * t26 + (t120 * t62 - t123 * t64) * qJD(4)) * t86, t149 * t123 - t71 * t168 + t230 * t25 + t64 * t84, -t149 * t120 - t71 * t167 + t230 * t26 - t62 * t84, -t230 * t48 + t71 * t84, (-t173 * t60 + t42) * t71 + t44 * t48 - (-t173 * t50 + t19) * t230 + t14 * t84 + t34 * t62 + t59 * t26 + t49 * t167 + ((-qJD(4) * t53 - t33) * t71 - t60 * t48 - (-qJD(4) * t32 - t20) * t230 + t150) * t120 + t155, -(-t174 * t60 + t169) * t71 - t196 * t48 + t161 * t230 - t15 * t84 + t34 * t64 - t59 * t25 - t49 * t168 + t150 * t123 + t156, -t1 * t230 + t12 * t48 + t120 * t154 + t167 * t27 + t22 * t62 + t26 * t35 + t3 * t71 + t8 * t84 + t155, -t10 * t84 + t123 * t154 - t16 * t48 - t168 * t27 + t2 * t230 + t22 * t64 - t25 * t35 - t4 * t71 + t156, t12 * t25 - t16 * t26 - t3 * t64 - t4 * t62 + (-t10 * t120 - t123 * t8) * t83 + t225 + (-t1 * t123 - t120 * t2 + (-t10 * t123 + t120 * t8) * qJD(4)) * t86, t1 * t12 + t10 * t4 + t2 * t16 + t27 * t22 + t8 * t3 + t5 * t35 + (-g(1) * t164 - g(2) * t137) * t125 + (g(1) * t137 - g(2) * t164) * t122; 0, 0, 0, -t170, -t159, -qJ(2) * t159 - t140, 0, 0, 0, 0, 0, t82 * t234 + t148, t233 * t234 + t135, 0, 0, 0, 0, 0, t6, t7, t6, t7, (t25 + t202) * t123 + t227 + t198, t153 * t120 + t228 * t123 - t27 * t82 - t223; 0, 0, 0, 0, 0, 0, -t82 * t233, -t233 ^ 2 + t82 ^ 2, t135, -t148, qJDD(3), -t82 * t92 + t139 + t144 - t206 + t232, -t233 * t92 - t147 + t220, t64 * t160 - t190, (-t25 + t202) * t123 - t227 + t198, t71 * t160 + t189 - t199, t141 + t201, -t71 * t82, -pkin(3) * t26 - t14 * t82 - t41 * t71 - t56 * t62 + (-t157 - t206) * t123 + (t226 * t71 + t138) * t120 + t193, pkin(3) * t25 + t197 * t71 + t15 * t82 - t56 * t64 + t99 + t138 * t123 + (-t139 + t157) * t120, -t105 * t26 - t28 * t62 + t48 * t95 - t8 * t82 + t194 * t71 + t166 * t123 + (-t27 * t233 + (t27 + t215) * qJD(4)) * t120 + t193, t10 * t82 + t105 * t25 - t28 * t64 + t48 * t96 + t99 - t195 * t71 + t27 * t160 + (pkin(4) * qJD(4) * t64 - t139 + t5) * t120, -t228 * t120 + t153 * t123 - t194 * t64 - t195 * t62 + t25 * t95 + t26 * t96 - t220, -t2 * t96 + t1 * t95 - t5 * t105 - g(3) * t143 + t194 * t8 + (pkin(4) * t174 - t28) * t27 + t195 * t10 + t152 * (t105 * t107 + t108 * t118); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t62, -t61 + t218, -t25 + t203, t200 - t26, t48, t15 * t71 - t49 * t64 + t129, t14 * t71 + t49 * t62 + t133, 0.2e1 * t216 + t192 + t204 + (t165 - t27) * t64 + t129, -pkin(4) * t218 + t191 + t71 * t9 + (qJD(5) + t27) * t62 + t133, pkin(4) * t25 - t217 * t62, t217 * t10 + (-t27 * t64 + t1 + t221) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 + t200, -t25 - t203, -t61 - t218, t10 * t62 + t64 * t8 - t139 - t166;];
tau_reg = t11;
