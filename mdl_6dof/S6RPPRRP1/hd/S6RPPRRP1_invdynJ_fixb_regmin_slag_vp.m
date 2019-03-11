% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:46
% EndTime: 2019-03-09 01:58:51
% DurationCPUTime: 2.21s
% Computational Cost: add. (3552->328), mult. (7635->408), div. (0->0), fcn. (5584->14), ass. (0->177)
t140 = sin(qJ(4));
t143 = cos(qJ(4));
t134 = sin(pkin(9));
t110 = pkin(1) * t134 + qJ(3);
t102 = t110 * qJD(1);
t135 = cos(pkin(10));
t120 = t135 * qJD(2);
t133 = sin(pkin(10));
t70 = t120 + (-pkin(7) * qJD(1) - t102) * t133;
t197 = qJD(1) * t135;
t83 = t133 * qJD(2) + t135 * t102;
t71 = pkin(7) * t197 + t83;
t32 = t140 * t70 + t143 * t71;
t244 = t32 * qJD(4);
t136 = cos(pkin(9));
t228 = pkin(1) * t136;
t114 = -pkin(2) - t228;
t192 = qJDD(1) * t114;
t100 = qJDD(3) + t192;
t132 = qJ(1) + pkin(9);
t122 = sin(t132);
t124 = cos(t132);
t182 = -g(1) * t122 + g(2) * t124;
t243 = -t100 - t182;
t118 = t135 * qJDD(2);
t96 = qJD(1) * qJD(3) + t110 * qJDD(1);
t63 = t118 + (-pkin(7) * qJDD(1) - t96) * t133;
t190 = t135 * qJDD(1);
t77 = t133 * qJDD(2) + t135 * t96;
t64 = pkin(7) * t190 + t77;
t159 = -t140 * t64 + t143 * t63 - t244;
t12 = -qJDD(4) * pkin(4) - t159;
t131 = pkin(10) + qJ(4);
t121 = sin(t131);
t123 = cos(t131);
t173 = g(1) * t124 + g(2) * t122;
t149 = -g(3) * t123 + t173 * t121;
t109 = t143 * t197;
t199 = t133 * t140;
t184 = qJD(1) * t199;
t90 = t109 - t184;
t85 = qJD(5) - t90;
t242 = -qJD(5) * pkin(8) * t85 - t12 + t149;
t139 = sin(qJ(5));
t178 = t139 * t85;
t142 = cos(qJ(5));
t99 = t133 * t143 + t135 * t140;
t91 = t99 * qJD(1);
t74 = qJD(4) * t139 + t142 * t91;
t241 = t74 * t178;
t239 = -t140 * t71 + t143 * t70;
t238 = t182 * t121;
t200 = t124 * t142;
t203 = t122 * t139;
t78 = t123 * t203 + t200;
t201 = t124 * t139;
t202 = t122 * t142;
t80 = -t123 * t201 + t202;
t237 = -g(1) * t80 + g(2) * t78;
t162 = t140 * t63 + t143 * t64;
t11 = qJDD(4) * pkin(8) + qJD(4) * t239 + t162;
t29 = qJD(4) * pkin(8) + t32;
t113 = pkin(3) * t135 + pkin(2);
t101 = -t113 - t228;
t87 = t101 * qJD(1) + qJD(3);
t37 = -t90 * pkin(4) - t91 * pkin(8) + t87;
t15 = t139 * t37 + t142 * t29;
t191 = t133 * qJDD(1);
t185 = qJD(4) * t109 + t140 * t190 + t143 * t191;
t152 = qJD(4) * t184 - t185;
t165 = t140 * t191 - t143 * t190;
t93 = t99 * qJD(4);
t57 = qJD(1) * t93 + t165;
t84 = t101 * qJDD(1) + qJDD(3);
t20 = t57 * pkin(4) + pkin(8) * t152 + t84;
t19 = t142 * t20;
t194 = t142 * qJD(4);
t196 = qJD(5) * t139;
t25 = -qJD(5) * t194 - t139 * qJDD(4) + t142 * t152 + t91 * t196;
t54 = qJDD(5) + t57;
t1 = t54 * pkin(5) + t25 * qJ(6) - t15 * qJD(5) - t74 * qJD(6) - t139 * t11 + t19;
t72 = t139 * t91 - t194;
t8 = -qJ(6) * t72 + t15;
t236 = t8 * t85 + t1;
t223 = g(3) * t121;
t235 = t173 * t123 + t223;
t98 = -t143 * t135 + t199;
t92 = t98 * qJD(4);
t166 = t54 * t99 - t85 * t92;
t187 = t99 * t196;
t234 = -t142 * t166 + t85 * t187;
t233 = t74 ^ 2;
t14 = -t139 * t29 + t142 * t37;
t7 = -qJ(6) * t74 + t14;
t6 = pkin(5) * t85 + t7;
t232 = -t7 + t6;
t215 = qJ(6) + pkin(8);
t181 = qJD(5) * t215;
t205 = qJ(6) * t142;
t55 = pkin(4) * t91 - pkin(8) * t90;
t48 = t142 * t55;
t229 = -pkin(5) * t91 - t142 * t181 + t90 * t205 - t48 + (-qJD(6) + t239) * t139;
t227 = pkin(5) * t139;
t141 = sin(qJ(1));
t221 = t141 * pkin(1);
t220 = t72 * t90;
t219 = t72 * t91;
t218 = t74 * t91;
t217 = t74 * t92;
t216 = pkin(7) + t110;
t148 = -t142 * qJDD(4) - t139 * t152;
t26 = t74 * qJD(5) + t148;
t214 = (-t26 * t99 + t72 * t92) * t142;
t195 = qJD(5) * t142;
t213 = -t139 * t26 - t72 * t195;
t212 = t139 * t55 + t142 * t239;
t211 = -t25 * t98 + t74 * t93;
t94 = t216 * t133;
t95 = t216 * t135;
t51 = -t140 * t94 + t143 * t95;
t43 = t142 * t51;
t44 = pkin(4) * t98 - pkin(8) * t99 + t101;
t210 = t139 * t44 + t43;
t206 = qJ(6) * t139;
t209 = t142 * qJD(6) - t139 * t181 + t90 * t206 - t212;
t208 = t139 * t54;
t204 = qJD(5) * t72;
t198 = t133 ^ 2 + t135 ^ 2;
t50 = t140 * t95 + t143 * t94;
t34 = -t98 * qJD(3) - t50 * qJD(4);
t56 = pkin(4) * t93 + pkin(8) * t92;
t189 = t139 * t56 + t142 * t34 + t44 * t195;
t186 = t99 * t195;
t183 = pkin(7) + qJ(3) + t227;
t180 = -qJD(5) * t37 - t11;
t179 = t142 * t85;
t176 = t74 * t186;
t175 = -t29 * t195 + t19;
t153 = t142 * t11 + t139 * t20 + t37 * t195 - t29 * t196;
t2 = -qJ(6) * t26 - qJD(6) * t72 + t153;
t174 = -t85 * t6 + t2;
t144 = cos(qJ(1));
t171 = g(1) * t141 - g(2) * t144;
t28 = -qJD(4) * pkin(4) - t239;
t170 = t12 * t99 - t28 * t92;
t169 = -t139 * t8 - t142 * t6;
t168 = t139 * t6 - t142 * t8;
t167 = -t98 * t26 - t93 * t72;
t76 = -t133 * t96 + t118;
t164 = -t76 * t133 + t77 * t135;
t163 = t133 * (-t102 * t133 + t120) - t135 * t83;
t160 = qJ(6) * t92 - qJD(6) * t99;
t116 = pkin(5) * t142 + pkin(4);
t158 = t116 * t123 + t121 * t215;
t156 = t142 * t54 + t90 * t178 - t85 * t196;
t155 = -pkin(8) * t54 + t85 * t28;
t154 = t113 + t158;
t151 = -t192 + t243;
t5 = t26 * pkin(5) + qJDD(6) + t12;
t147 = -t166 * t139 - t85 * t186;
t35 = t99 * qJD(3) + t51 * qJD(4);
t128 = t144 * pkin(1);
t104 = t215 * t142;
t103 = t215 * t139;
t81 = t123 * t200 + t203;
t79 = -t123 * t202 + t201;
t69 = t72 ^ 2;
t59 = -qJD(4) * t93 - qJDD(4) * t98;
t58 = -qJD(4) * t92 + qJDD(4) * t99;
t49 = t142 * t56;
t41 = t142 * t44;
t21 = pkin(5) * t72 + qJD(6) + t28;
t17 = -t99 * t206 + t210;
t16 = pkin(5) * t98 - t139 * t51 - t99 * t205 + t41;
t4 = -qJ(6) * t186 + (-qJD(5) * t51 + t160) * t139 + t189;
t3 = t93 * pkin(5) - t139 * t34 + t49 + t160 * t142 + (-t43 + (qJ(6) * t99 - t44) * t139) * qJD(5);
t9 = [qJDD(1), t171, g(1) * t144 + g(2) * t141 (t171 + (t134 ^ 2 + t136 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t151 * t135, -t151 * t133, t96 * t198 + t164 - t173, t100 * t114 - g(1) * (-pkin(2) * t122 + qJ(3) * t124 - t221) - g(2) * (pkin(2) * t124 + qJ(3) * t122 + t128) + t164 * t110 - t163 * qJD(3), -t152 * t99 - t91 * t92, t152 * t98 - t99 * t57 - t92 * t90 - t91 * t93, t58, t59, 0, -t35 * qJD(4) - t50 * qJDD(4) + t101 * t57 - t123 * t182 + t84 * t98 + t87 * t93, -t34 * qJD(4) - t51 * qJDD(4) - t101 * t152 + t84 * t99 - t87 * t92 + t238, -t74 * t187 + (-t25 * t99 - t217) * t142, -t176 + (t217 + (t25 + t204) * t99) * t139 + t214, t211 - t234, t147 + t167, t54 * t98 + t85 * t93 (-t195 * t51 + t49) * t85 + t41 * t54 + t175 * t98 + t14 * t93 + t35 * t72 + t50 * t26 + t28 * t186 - g(1) * t79 - g(2) * t81 + ((-qJD(5) * t44 - t34) * t85 - t51 * t54 + t180 * t98 + t170) * t139 -(-t196 * t51 + t189) * t85 - t210 * t54 - t153 * t98 - t15 * t93 + t35 * t74 - t50 * t25 - t28 * t187 - g(1) * t78 - g(2) * t80 + t170 * t142, t16 * t25 - t17 * t26 - t3 * t74 - t4 * t72 - t169 * t92 - t238 + (qJD(5) * t168 - t1 * t142 - t2 * t139) * t99, t2 * t17 + t8 * t4 + t1 * t16 + t6 * t3 + t5 * (t99 * t227 + t50) + t21 * ((-t139 * t92 + t186) * pkin(5) + t35) + g(1) * t221 - g(2) * t128 + (-g(1) * t183 - g(2) * t154) * t124 + (g(1) * t154 - g(2) * t183) * t122; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, t133 * t77 + t135 * t76 - g(3), 0, 0, 0, 0, 0, t59, -t58, 0, 0, 0, 0, 0, t147 - t167, t211 + t234, t176 + (-t217 + (-t25 + t204) * t99) * t139 + t214, t21 * t93 + t5 * t98 - g(3) + t168 * t92 + (qJD(5) * t169 - t1 * t139 + t2 * t142) * t99; 0, 0, 0, 0, -t190, t191, -t198 * qJD(1) ^ 2, qJD(1) * t163 - t243, 0, 0, 0, 0, 0, 0.2e1 * t91 * qJD(4) + t165 (t90 - t184) * qJD(4) + t185, 0, 0, 0, 0, 0, t156 - t219, -t142 * t85 ^ 2 - t208 - t218 (t25 + t220) * t142 + t241 + t213, t174 * t139 + t236 * t142 - t21 * t91 + t182; 0, 0, 0, 0, 0, 0, 0, 0, -t91 * t90, -t90 ^ 2 + t91 ^ 2 (-t90 - t184) * qJD(4) + t185, -t165, qJDD(4), -t87 * t91 + t149 + t159 + t244, -t87 * t90 - t162 + t235, -t25 * t139 + t179 * t74 (-t25 + t220) * t142 - t241 + t213, t179 * t85 + t208 - t218, t156 + t219, -t85 * t91, -pkin(4) * t26 - t14 * t91 - t32 * t72 - t48 * t85 + (t239 * t85 + t155) * t139 + t242 * t142, pkin(4) * t25 - t242 * t139 + t155 * t142 + t15 * t91 + t212 * t85 - t32 * t74, -t103 * t25 - t104 * t26 - t236 * t139 + t174 * t142 - t209 * t72 - t229 * t74 - t235, t2 * t104 - t1 * t103 - t5 * t116 - g(3) * t158 + t209 * t8 + t229 * t6 + (pkin(5) * t178 - t32) * t21 + t173 * (t116 * t121 - t123 * t215); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t72, -t69 + t233, t72 * t85 - t25, -t148 + (-qJD(5) + t85) * t74, t54, t15 * t85 - t28 * t74 + (t180 + t223) * t139 + t175 + t237, g(1) * t81 - g(2) * t79 + t14 * t85 + t142 * t223 + t28 * t72 - t153, pkin(5) * t25 - t232 * t72, t232 * t8 + (t139 * t223 - t21 * t74 + t1 + t237) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 - t233, t6 * t74 + t8 * t72 - t149 + t5;];
tau_reg  = t9;
