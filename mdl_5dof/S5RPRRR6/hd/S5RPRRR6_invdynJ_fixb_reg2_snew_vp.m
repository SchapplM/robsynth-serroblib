% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:41
% EndTime: 2019-12-31 19:01:49
% DurationCPUTime: 2.52s
% Computational Cost: add. (8548->315), mult. (17019->441), div. (0->0), fcn. (11432->10), ass. (0->191)
t171 = sin(qJ(5));
t172 = sin(qJ(4));
t173 = sin(qJ(3));
t176 = cos(qJ(4));
t177 = cos(qJ(3));
t136 = (t177 * t172 + t173 * t176) * qJD(1);
t163 = qJD(3) + qJD(4);
t175 = cos(qJ(5));
t117 = t171 * t136 - t175 * t163;
t119 = t175 * t136 + t171 * t163;
t96 = t119 * t117;
t207 = qJD(1) * qJD(3);
t197 = t177 * t207;
t206 = t173 * qJDD(1);
t141 = t197 + t206;
t159 = t177 * qJDD(1);
t198 = t173 * t207;
t142 = t159 - t198;
t192 = t172 * t141 - t176 * t142;
t101 = -t136 * qJD(4) - t192;
t98 = qJDD(5) - t101;
t233 = -t96 + t98;
t237 = t171 * t233;
t210 = qJD(1) * t173;
t134 = -t176 * t177 * qJD(1) + t172 * t210;
t114 = t136 * t134;
t162 = qJDD(3) + qJDD(4);
t232 = -t114 + t162;
t236 = t172 * t232;
t235 = t175 * t233;
t234 = t176 * t232;
t112 = t134 * pkin(4) - t136 * pkin(8);
t230 = t163 ^ 2;
t166 = -g(3) + qJDD(2);
t180 = qJD(1) ^ 2;
t174 = sin(qJ(1));
t178 = cos(qJ(1));
t191 = t178 * g(1) + t174 * g(2);
t139 = -t180 * pkin(1) - t191;
t168 = sin(pkin(9));
t169 = cos(pkin(9));
t190 = t174 * g(1) - t178 * g(2);
t185 = qJDD(1) * pkin(1) + t190;
t211 = t169 * t139 + t168 * t185;
t186 = -t180 * pkin(2) + qJDD(1) * pkin(6) + t211;
t100 = t173 * t166 + t177 * t186;
t148 = qJD(3) * pkin(3) - pkin(7) * t210;
t165 = t177 ^ 2;
t161 = t165 * t180;
t81 = -pkin(3) * t161 + t142 * pkin(7) - qJD(3) * t148 + t100;
t222 = t176 * t81;
t183 = t173 * t186;
t214 = t173 * t180;
t227 = t141 * pkin(7);
t231 = qJDD(3) * pkin(3) + (pkin(3) * t214 + pkin(7) * t207 + t166) * t177 - t183 - t227;
t49 = t231 * t172 + t222;
t38 = -t230 * pkin(4) + t162 * pkin(8) - t134 * t112 + t49;
t193 = -t168 * t139 + t169 * t185;
t107 = -qJDD(1) * pkin(2) - t180 * pkin(6) - t193;
t84 = -t142 * pkin(3) - pkin(7) * t161 + t148 * t210 + t107;
t102 = -t134 * qJD(4) + t176 * t141 + t172 * t142;
t126 = t163 * t134;
t90 = t102 - t126;
t40 = -t90 * pkin(8) + (t163 * t136 - t101) * pkin(4) + t84;
t16 = t171 * t38 - t175 * t40;
t17 = t171 * t40 + t175 * t38;
t7 = t171 * t16 + t175 * t17;
t128 = qJD(5) + t134;
t194 = t171 * t102 - t175 * t162;
t64 = (qJD(5) - t128) * t119 + t194;
t115 = t117 ^ 2;
t116 = t119 ^ 2;
t127 = t128 ^ 2;
t129 = t134 ^ 2;
t130 = t136 ^ 2;
t48 = t172 * t81 - t176 * t231;
t37 = -t162 * pkin(4) - t230 * pkin(8) + t136 * t112 + t48;
t229 = -pkin(4) * t37 + pkin(8) * t7;
t228 = pkin(4) * t172;
t34 = t171 * t37;
t72 = t96 + t98;
t226 = t171 * t72;
t225 = t172 * t84;
t22 = t172 * t49 - t176 * t48;
t224 = t173 * t22;
t35 = t175 * t37;
t223 = t175 * t72;
t221 = t176 * t84;
t220 = t128 * t171;
t219 = t128 * t175;
t218 = t163 * t172;
t217 = t163 * t176;
t110 = t114 + t162;
t216 = t172 * t110;
t152 = t177 * t214;
t146 = qJDD(3) + t152;
t215 = t173 * t146;
t213 = t176 * t110;
t147 = qJDD(3) - t152;
t212 = t177 * t147;
t208 = qJD(5) + t128;
t94 = -t116 - t127;
t46 = -t171 * t94 - t223;
t187 = -t175 * t102 - t171 * t162;
t69 = t208 * t117 + t187;
t205 = pkin(4) * t69 + pkin(8) * t46 + t34;
t85 = -t127 - t115;
t43 = t175 * t85 - t237;
t65 = -t208 * t119 - t194;
t204 = pkin(4) * t65 + pkin(8) * t43 - t35;
t203 = t172 * t96;
t202 = t176 * t96;
t201 = -pkin(1) * t169 - pkin(2);
t200 = pkin(1) * t168 + pkin(6);
t199 = -pkin(4) * t176 - pkin(3);
t23 = t172 * t48 + t176 * t49;
t99 = -t177 * t166 + t183;
t70 = t177 * t100 + t173 * t99;
t106 = t128 * t117;
t77 = -t117 * qJD(5) - t187;
t68 = t106 + t77;
t33 = t171 * t68 - t175 * t64;
t80 = t115 + t116;
t196 = pkin(4) * t80 + pkin(8) * t33 + t7;
t188 = -t175 * t16 + t171 * t17;
t184 = (-qJD(4) + t163) * t136 - t192;
t179 = qJD(3) ^ 2;
t164 = t173 ^ 2;
t160 = t164 * t180;
t150 = -t161 - t179;
t149 = -t160 - t179;
t145 = t160 + t161;
t144 = (t164 + t165) * qJDD(1);
t143 = t159 - 0.2e1 * t198;
t140 = 0.2e1 * t197 + t206;
t124 = -t130 + t230;
t123 = t129 - t230;
t122 = -t130 - t230;
t121 = -t173 * t149 - t212;
t120 = t177 * t150 - t215;
t113 = t130 - t129;
t108 = -t230 - t129;
t105 = -t116 + t127;
t104 = t115 - t127;
t103 = -t129 - t130;
t95 = t116 - t115;
t93 = -t172 * t122 - t213;
t92 = t176 * t122 - t216;
t91 = t102 + t126;
t86 = (qJD(4) + t163) * t136 + t192;
t83 = t176 * t108 - t236;
t82 = t172 * t108 + t234;
t76 = -t119 * qJD(5) - t194;
t75 = (-t117 * t175 + t119 * t171) * t128;
t74 = (-t117 * t171 - t119 * t175) * t128;
t67 = -t106 + t77;
t61 = -t119 * t220 + t175 * t77;
t60 = t119 * t219 + t171 * t77;
t59 = t117 * t219 - t171 * t76;
t58 = t117 * t220 + t175 * t76;
t57 = -t173 * t92 + t177 * t93;
t56 = t172 * t91 + t176 * t184;
t55 = t172 * t184 - t176 * t91;
t54 = t175 * t104 - t226;
t53 = -t171 * t105 + t235;
t52 = t171 * t104 + t223;
t51 = t175 * t105 + t237;
t50 = -t173 * t82 + t177 * t83;
t45 = t175 * t94 - t226;
t42 = t171 * t85 + t235;
t32 = -t171 * t67 + t175 * t65;
t31 = -t171 * t64 - t175 * t68;
t30 = t171 * t65 + t175 * t67;
t28 = -t173 * t55 + t177 * t56;
t27 = -t172 * t69 + t176 * t46;
t26 = t172 * t46 + t176 * t69;
t25 = -t172 * t65 + t176 * t43;
t24 = t172 * t43 + t176 * t65;
t21 = -t172 * t80 + t176 * t33;
t20 = t172 * t33 + t176 * t80;
t19 = -pkin(8) * t45 + t35;
t18 = -pkin(8) * t42 + t34;
t13 = -pkin(4) * t45 + t17;
t12 = -t173 * t26 + t177 * t27;
t11 = -pkin(4) * t42 + t16;
t10 = -t173 * t24 + t177 * t25;
t9 = t177 * t23 - t224;
t4 = t172 * t37 + t176 * t7;
t3 = t172 * t7 - t176 * t37;
t2 = -pkin(8) * t31 - t188;
t1 = [0, 0, 0, 0, 0, qJDD(1), t190, t191, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t169 * qJDD(1) - t168 * t180) + t193, pkin(1) * (-t168 * qJDD(1) - t169 * t180) - t211, 0, pkin(1) * (t168 * t211 + t169 * t193), (t141 + t197) * t173, t177 * t140 + t173 * t143, t215 + t177 * (-t160 + t179), (t142 - t198) * t177, t173 * (t161 - t179) + t212, 0, -t177 * t107 + pkin(2) * t143 + pkin(6) * t120 + pkin(1) * (t168 * t120 + t169 * t143), t173 * t107 - pkin(2) * t140 + pkin(6) * t121 + pkin(1) * (t168 * t121 - t169 * t140), pkin(2) * t145 + pkin(6) * t144 + pkin(1) * (t168 * t144 + t169 * t145) + t70, -pkin(2) * t107 + pkin(6) * t70 + pkin(1) * (-t169 * t107 + t168 * t70), t173 * (t176 * t102 - t136 * t218) + t177 * (t172 * t102 + t136 * t217), t173 * (-t172 * t90 - t176 * t86) + t177 * (-t172 * t86 + t176 * t90), t173 * (-t172 * t124 + t234) + t177 * (t176 * t124 + t236), t173 * (-t172 * t101 + t134 * t217) + t177 * (t176 * t101 + t134 * t218), t173 * (t176 * t123 - t216) + t177 * (t172 * t123 + t213), (t173 * (-t134 * t176 + t136 * t172) + t177 * (-t134 * t172 - t136 * t176)) * t163, t173 * (-pkin(7) * t82 + t225) + t177 * (-pkin(3) * t86 + pkin(7) * t83 - t221) - pkin(2) * t86 + pkin(6) * t50 + pkin(1) * (t168 * t50 - t169 * t86), t173 * (-pkin(7) * t92 + t221) + t177 * (-pkin(3) * t90 + pkin(7) * t93 + t225) - pkin(2) * t90 + pkin(6) * t57 + pkin(1) * (t168 * t57 - t169 * t90), t173 * (-pkin(7) * t55 - t22) + t177 * (-pkin(3) * t103 + pkin(7) * t56 + t23) - pkin(2) * t103 + pkin(6) * t28 + pkin(1) * (-t169 * t103 + t168 * t28), -pkin(7) * t224 + t177 * (-pkin(3) * t84 + pkin(7) * t23) - pkin(2) * t84 + pkin(6) * t9 + pkin(1) * (t168 * t9 - t169 * t84), t173 * (t176 * t61 + t203) + t177 * (t172 * t61 - t202), t173 * (t172 * t95 + t176 * t32) + t177 * (t172 * t32 - t176 * t95), t173 * (t172 * t68 + t176 * t53) + t177 * (t172 * t53 - t176 * t68), t173 * (t176 * t59 - t203) + t177 * (t172 * t59 + t202), t173 * (-t172 * t64 + t176 * t54) + t177 * (t172 * t54 + t176 * t64), t173 * (t172 * t98 + t176 * t75) + t177 * (t172 * t75 - t176 * t98), t173 * (-pkin(7) * t24 - t172 * t11 + t176 * t18) + t177 * (-pkin(3) * t42 + pkin(7) * t25 + t176 * t11 + t172 * t18) - pkin(2) * t42 + pkin(6) * t10 + pkin(1) * (t168 * t10 - t169 * t42), t173 * (-pkin(7) * t26 - t172 * t13 + t176 * t19) + t177 * (-pkin(3) * t45 + pkin(7) * t27 + t176 * t13 + t172 * t19) - pkin(2) * t45 + pkin(6) * t12 + pkin(1) * (t168 * t12 - t169 * t45), t173 * (-pkin(7) * t20 + t176 * t2) + t177 * (pkin(7) * t21 + t172 * t2) + t200 * (-t173 * t20 + t177 * t21) + (t173 * t228 + t177 * t199 + t201) * t31, (t173 * (-pkin(8) * t176 + t228) + t177 * (-pkin(8) * t172 + t199) + t201) * t188 + (t200 + pkin(7)) * (-t173 * t3 + t177 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, 0, 0, 0, 0, 0, 0, t177 * t146 + t173 * t150, -t173 * t147 + t177 * t149, 0, t173 * t100 - t177 * t99, 0, 0, 0, 0, 0, 0, t173 * t83 + t177 * t82, t173 * t93 + t177 * t92, t173 * t56 + t177 * t55, t173 * t23 + t177 * t22, 0, 0, 0, 0, 0, 0, t173 * t25 + t177 * t24, t173 * t27 + t177 * t26, t173 * t21 + t177 * t20, t173 * t4 + t177 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, t160 - t161, t206, t152, t159, qJDD(3), -t99, -t100, 0, 0, t114, t113, t91, -t114, t184, t162, pkin(3) * t82 - t48, -t222 - t172 * (pkin(7) * t197 - t227 - t99) + (-t172 * t146 + t92) * pkin(3), pkin(3) * t55, pkin(3) * t22, t60, t30, t51, t58, t52, t74, pkin(3) * t24 + t204, pkin(3) * t26 + t205, pkin(3) * t20 + t196, pkin(3) * t3 + t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t113, t91, -t114, t184, t162, -t48, -t49, 0, 0, t60, t30, t51, t58, t52, t74, t204, t205, t196, t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t95, t68, -t96, -t64, t98, -t16, -t17, 0, 0;];
tauJ_reg = t1;
