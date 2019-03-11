% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [6x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:06
% EndTime: 2019-03-09 03:03:12
% DurationCPUTime: 2.18s
% Computational Cost: add. (3907->351), mult. (8426->451), div. (0->0), fcn. (5856->14), ass. (0->194)
t139 = sin(pkin(9));
t117 = t139 * pkin(1) + pkin(7);
t212 = qJ(4) + t117;
t138 = sin(pkin(10));
t140 = cos(pkin(10));
t145 = sin(qJ(3));
t148 = cos(qJ(3));
t102 = t138 * t148 + t140 * t145;
t141 = cos(pkin(9));
t119 = -t141 * pkin(1) - pkin(2);
t132 = t148 * pkin(3);
t257 = t119 - t132;
t259 = -t102 * pkin(8) + t257;
t144 = sin(qJ(5));
t207 = qJD(1) * t145;
t213 = t140 * t148;
t92 = qJD(1) * t213 - t138 * t207;
t90 = qJD(5) - t92;
t186 = t144 * t90;
t147 = cos(qJ(5));
t94 = t102 * qJD(1);
t77 = t144 * qJD(3) + t147 * t94;
t258 = t77 * t186;
t206 = qJD(5) * t144;
t201 = t148 * qJDD(1);
t202 = t145 * qJDD(1);
t173 = -t138 * t202 + t140 * t201;
t93 = t102 * qJD(3);
t65 = qJD(1) * t93 + qJDD(5) - t173;
t166 = -t147 * t65 + t90 * t206;
t135 = qJ(1) + pkin(9);
t124 = sin(t135);
t126 = cos(t135);
t189 = -g(1) * t124 + g(2) * t126;
t179 = g(1) * t126 + g(2) * t124;
t134 = qJ(3) + pkin(10);
t125 = cos(t134);
t214 = t126 * t147;
t217 = t124 * t144;
t85 = t125 * t217 + t214;
t215 = t126 * t144;
t216 = t124 * t147;
t87 = -t125 * t215 + t216;
t256 = -g(1) * t87 + g(2) * t85;
t129 = t148 * qJDD(2);
t106 = t117 * qJDD(1);
t159 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(3) * qJD(2) + t106;
t183 = t212 * qJD(1);
t168 = t183 * qJD(3);
t34 = qJDD(3) * pkin(3) - t159 * t145 - t148 * t168 + t129;
t38 = (qJDD(2) - t168) * t145 + t159 * t148;
t13 = -t138 * t38 + t140 * t34;
t11 = -qJDD(3) * pkin(4) - t13;
t116 = t138 * pkin(3) + pkin(8);
t123 = sin(t134);
t160 = g(3) * t125 - t179 * t123;
t255 = qJD(5) * t116 * t90 + t11 + t160;
t14 = t138 * t34 + t140 * t38;
t12 = qJDD(3) * pkin(8) + t14;
t80 = t145 * qJD(2) + t183 * t148;
t227 = t140 * t80;
t228 = qJD(3) * pkin(3);
t79 = t148 * qJD(2) - t183 * t145;
t72 = t79 + t228;
t37 = t138 * t72 + t227;
t31 = qJD(3) * pkin(8) + t37;
t91 = qJD(1) * t257 + qJD(4);
t46 = -t92 * pkin(4) - t94 * pkin(8) + t91;
t17 = t144 * t46 + t147 * t31;
t203 = qJD(1) * qJD(3);
t192 = t145 * t203;
t110 = t138 * t192;
t191 = t148 * t203;
t167 = t140 * t191 - t110;
t200 = pkin(3) * t192 + qJDD(4);
t66 = -qJD(3) * t94 + t173;
t22 = -t66 * pkin(4) - t167 * pkin(8) + qJDD(1) * t259 + t200;
t21 = t147 * t22;
t204 = t147 * qJD(3);
t252 = t102 * qJDD(1) + t167;
t27 = -qJD(5) * t204 - t144 * qJDD(3) - t147 * t252 + t94 * t206;
t1 = t65 * pkin(5) + t27 * qJ(6) - t17 * qJD(5) - t77 * qJD(6) - t144 * t12 + t21;
t75 = t144 * t94 - t204;
t8 = -t75 * qJ(6) + t17;
t254 = t8 * t90 + t1;
t101 = t138 * t145 - t213;
t96 = t101 * qJD(3);
t221 = t147 * t96;
t253 = -t166 * t102 - t90 * t221;
t251 = t77 ^ 2;
t16 = -t144 * t31 + t147 * t46;
t7 = -t77 * qJ(6) + t16;
t6 = t90 * pkin(5) + t7;
t250 = -t7 + t6;
t211 = qJ(6) + t116;
t184 = qJD(5) * t211;
t69 = t138 * t80;
t40 = t140 * t79 - t69;
t55 = pkin(3) * t207 + t94 * pkin(4) - t92 * pkin(8);
t49 = t147 * t55;
t247 = -t94 * pkin(5) - t49 + (qJ(6) * t92 - t184) * t147 + (-qJD(6) + t40) * t144;
t243 = g(3) * t123;
t241 = g(3) * t148;
t239 = t140 * pkin(3);
t146 = sin(qJ(1));
t238 = t146 * pkin(1);
t237 = t75 * t92;
t236 = t77 * t94;
t235 = t94 * t75;
t218 = t102 * t147;
t153 = -t147 * qJDD(3) + t144 * t252;
t28 = t77 * qJD(5) + t153;
t234 = -t28 * t218 + t75 * t221;
t205 = qJD(5) * t147;
t233 = -t144 * t28 - t75 * t205;
t232 = t144 * t55 + t147 * t40;
t231 = -t27 * t101 + t77 * t93;
t100 = t212 * t148;
t98 = t212 * t145;
t63 = t140 * t100 - t138 * t98;
t58 = t147 * t63;
t59 = t101 * pkin(4) + t259;
t230 = t144 * t59 + t58;
t224 = t144 * t92;
t229 = qJ(6) * t224 + t147 * qJD(6) - t144 * t184 - t232;
t226 = t144 * t65;
t225 = t144 * t75;
t223 = t144 * t96;
t222 = t147 * t77;
t220 = t27 * t144;
t219 = t102 * t144;
t210 = qJDD(2) - g(3);
t122 = t132 + pkin(2);
t149 = cos(qJ(1));
t209 = t149 * pkin(1) + t126 * t122;
t136 = t145 ^ 2;
t208 = -t148 ^ 2 + t136;
t109 = qJD(1) * t119;
t185 = qJD(3) * t212;
t81 = t148 * qJD(4) - t145 * t185;
t82 = -t145 * qJD(4) - t148 * t185;
t45 = t138 * t82 + t140 * t81;
t196 = t145 * t228;
t56 = t93 * pkin(4) + t96 * pkin(8) + t196;
t199 = t144 * t56 + t147 * t45 + t59 * t205;
t197 = t77 * t223;
t121 = t147 * pkin(5) + pkin(4);
t193 = t102 * t205;
t143 = -qJ(4) - pkin(7);
t190 = pkin(5) * t144 - t143;
t39 = t138 * t79 + t227;
t44 = t138 * t81 - t140 * t82;
t36 = t140 * t72 - t69;
t188 = -qJD(5) * t46 - t12;
t62 = t138 * t100 + t140 * t98;
t187 = t147 * t90;
t181 = -t31 * t205 + t21;
t165 = t147 * t12 + t144 * t22 + t46 * t205 - t31 * t206;
t2 = -t28 * qJ(6) - t75 * qJD(6) + t165;
t180 = -t90 * t6 + t2;
t177 = g(1) * t146 - g(2) * t149;
t176 = -t144 * t8 - t147 * t6;
t175 = t144 * t6 - t147 * t8;
t174 = -t101 * t28 - t93 * t75;
t172 = qJ(6) * t96 - qJD(6) * t102;
t142 = -qJ(6) - pkin(8);
t171 = t125 * t121 - t123 * t142;
t169 = t90 * t224 - t166;
t30 = -qJD(3) * pkin(4) - t36;
t164 = -t116 * t65 + t90 * t30;
t162 = -t109 * qJD(1) - t106 + t179;
t161 = 0.2e1 * t109 * qJD(3) - qJDD(3) * t117;
t5 = t28 * pkin(5) + qJDD(6) + t11;
t158 = qJDD(1) * t257 + t200;
t150 = qJD(3) ^ 2;
t157 = -0.2e1 * qJDD(1) * t119 - t117 * t150 - t189;
t156 = (-t90 * t205 - t226) * t102 + t90 * t223;
t151 = qJD(1) ^ 2;
t118 = -pkin(4) - t239;
t105 = qJDD(3) * t148 - t150 * t145;
t104 = qJDD(3) * t145 + t150 * t148;
t99 = t211 * t147;
t97 = t211 * t144;
t88 = t125 * t214 + t217;
t86 = -t125 * t216 + t215;
t74 = t75 ^ 2;
t54 = t147 * t59;
t50 = t147 * t56;
t23 = t75 * pkin(5) + qJD(6) + t30;
t19 = -qJ(6) * t219 + t230;
t18 = t101 * pkin(5) - qJ(6) * t218 - t144 * t63 + t54;
t4 = -qJ(6) * t193 + (-qJD(5) * t63 + t172) * t144 + t199;
t3 = t93 * pkin(5) - t144 * t45 + t50 + t172 * t147 + (-t58 + (qJ(6) * t102 - t59) * t144) * qJD(5);
t9 = [qJDD(1), t177, g(1) * t149 + g(2) * t146 (t177 + (t139 ^ 2 + t141 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t136 * qJDD(1) + 0.2e1 * t145 * t191, 0.2e1 * t145 * t201 - 0.2e1 * t208 * t203, t104, t105, 0, t145 * t161 + t148 * t157, -t145 * t157 + t148 * t161, -t14 * t101 - t13 * t102 + t252 * t62 + t36 * t96 - t37 * t93 + t44 * t94 + t45 * t92 + t63 * t66 - t179, t14 * t63 + t37 * t45 - t13 * t62 - t36 * t44 + t158 * t257 + t91 * t196 - g(1) * (-t124 * t122 - t126 * t143 - t238) - g(2) * (-t124 * t143 + t209) -t77 * t221 + (-t147 * t27 - t77 * t206) * t102, t197 + (t220 + (-t222 + t225) * qJD(5)) * t102 + t234, t231 + t253, t156 + t174, t65 * t101 + t90 * t93 (-t63 * t205 + t50) * t90 + t54 * t65 + t181 * t101 + t16 * t93 + t44 * t75 + t62 * t28 + t30 * t193 - g(1) * t86 - g(2) * t88 + ((-qJD(5) * t59 - t45) * t90 - t63 * t65 + t188 * t101 + t11 * t102 - t30 * t96) * t144 -(-t206 * t63 + t199) * t90 - t230 * t65 - t165 * t101 - t17 * t93 + t44 * t77 - t62 * t27 - t30 * t221 - g(1) * t85 - g(2) * t87 + (t11 * t147 - t206 * t30) * t102, t18 * t27 - t19 * t28 - t3 * t77 - t4 * t75 - t176 * t96 - t189 * t123 + (qJD(5) * t175 - t1 * t147 - t2 * t144) * t102, t2 * t19 + t8 * t4 + t1 * t18 + t6 * t3 + t5 * (pkin(5) * t219 + t62) + t23 * ((t193 - t223) * pkin(5) + t44) + g(1) * t238 - g(2) * t209 + (-g(1) * t190 - g(2) * t171) * t126 + (-g(1) * (-t122 - t171) - g(2) * t190) * t124; 0, 0, 0, t210, 0, 0, 0, 0, 0, t105, -t104, t101 * t252 + t102 * t66 - t96 * t92 + t93 * t94, -t13 * t101 + t14 * t102 - t36 * t93 - t37 * t96 - g(3), 0, 0, 0, 0, 0, t156 - t174, t231 - t253, -t197 + (-t220 + (t222 + t225) * qJD(5)) * t102 + t234, t5 * t101 + t23 * t93 - g(3) + t175 * t96 + (qJD(5) * t176 - t1 * t144 + t2 * t147) * t102; 0, 0, 0, 0, -t145 * t151 * t148, t208 * t151, t202, t201, qJDD(3), t145 * t162 + t129 - t241, -t210 * t145 + t162 * t148 (t37 - t39) * t94 + (-t40 + t36) * t92 + (t138 * t66 + (-t138 * t201 + t110 + (-t191 - t202) * t140) * t140) * pkin(3), t36 * t39 - t37 * t40 + (-t241 + t13 * t140 + t138 * t14 + (-qJD(1) * t91 + t179) * t145) * pkin(3), t77 * t187 - t220 (-t27 + t237) * t147 - t258 + t233, t90 * t187 + t226 - t236, t169 + t235, -t90 * t94, t118 * t28 - t16 * t94 - t39 * t75 - t49 * t90 + (t40 * t90 + t164) * t144 - t255 * t147, -t118 * t27 + t255 * t144 + t164 * t147 + t17 * t94 + t232 * t90 - t39 * t77, -t179 * t125 - t254 * t144 + t180 * t147 - t229 * t75 - t247 * t77 - t97 * t27 - t99 * t28 - t243, t2 * t99 - t1 * t97 + t5 * (-t121 - t239) - g(3) * (t132 + t171) + t229 * t8 + t247 * t6 + (pkin(5) * t186 - t39) * t23 + t179 * (pkin(3) * t145 + t121 * t123 + t125 * t142); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92 ^ 2 - t94 ^ 2, t36 * t94 - t37 * t92 + t158 + t189, 0, 0, 0, 0, 0, t169 - t235, -t147 * t90 ^ 2 - t226 - t236 (t27 + t237) * t147 + t258 + t233, t180 * t144 + t254 * t147 - t23 * t94 + t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t75, -t74 + t251, t75 * t90 - t27, -t153 + (-qJD(5) + t90) * t77, t65, t17 * t90 - t30 * t77 + (t188 + t243) * t144 + t181 + t256, g(1) * t88 - g(2) * t86 + t147 * t243 + t16 * t90 + t30 * t75 - t165, pkin(5) * t27 - t250 * t75, t250 * t8 + (t144 * t243 - t23 * t77 + t1 + t256) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74 - t251, t6 * t77 + t8 * t75 + t160 + t5;];
tau_reg  = t9;
