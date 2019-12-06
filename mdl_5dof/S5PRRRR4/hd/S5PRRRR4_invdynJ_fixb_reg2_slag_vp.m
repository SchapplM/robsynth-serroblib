% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:08:04
% EndTime: 2019-12-05 17:08:07
% DurationCPUTime: 1.61s
% Computational Cost: add. (2592->266), mult. (3970->339), div. (0->0), fcn. (2480->12), ass. (0->178)
t126 = qJDD(2) + qJDD(3);
t136 = sin(qJ(3));
t191 = qJDD(2) * t136;
t139 = cos(qJ(3));
t198 = qJD(3) * t139;
t62 = t126 * pkin(7) + (qJD(2) * t198 + t191) * pkin(2);
t239 = qJD(1) * qJD(4) + t62;
t127 = pkin(9) + qJ(2);
t119 = qJ(3) + t127;
t107 = cos(t119);
t230 = g(2) * t107;
t199 = qJD(3) * t136;
t188 = pkin(2) * t199;
t233 = pkin(2) * t139;
t203 = -qJD(2) * t188 + qJDD(2) * t233;
t228 = t126 * pkin(3);
t61 = -t203 - t228;
t238 = t61 + t230;
t137 = cos(qJ(5));
t138 = cos(qJ(4));
t194 = qJD(5) * t137;
t196 = qJD(4) * t138;
t237 = -t137 * t196 - t138 * t194;
t129 = qJD(2) + qJD(3);
t134 = sin(qJ(5));
t135 = sin(qJ(4));
t206 = t135 * t137;
t70 = t134 * t138 + t206;
t59 = t70 * t129;
t106 = sin(t119);
t236 = g(1) * t107 + g(2) * t106;
t184 = t135 * qJDD(1) + t138 * t239;
t197 = qJD(4) * t135;
t220 = pkin(2) * qJD(2);
t189 = t136 * t220;
t81 = pkin(7) * t129 + t189;
t31 = -t197 * t81 + t184;
t118 = t138 * qJDD(1);
t193 = t135 * qJD(1);
t215 = t138 * t81;
t56 = t193 + t215;
t32 = -qJD(4) * t56 - t135 * t62 + t118;
t120 = t138 * qJD(1);
t218 = t135 * t81;
t55 = t120 - t218;
t146 = -t32 * t135 + t31 * t138 + (-t135 * t56 - t138 * t55) * qJD(4);
t130 = t135 ^ 2;
t131 = t138 ^ 2;
t201 = t130 + t131;
t235 = t129 * t201;
t128 = qJD(4) + qJD(5);
t176 = pkin(8) * t129 + t81;
t46 = t138 * t176 + t193;
t16 = qJDD(4) * pkin(4) + t118 + (-pkin(8) * t126 - t62) * t135 - t46 * qJD(4);
t179 = t129 * t197;
t204 = t138 * t126;
t154 = t179 - t204;
t17 = -pkin(8) * t154 + t31;
t195 = qJD(5) * t134;
t45 = -t135 * t176 + t120;
t41 = qJD(4) * pkin(4) + t45;
t4 = (qJD(5) * t41 + t17) * t137 + t134 * t16 - t46 * t195;
t159 = t135 * t55 - t138 * t56;
t200 = qJD(2) * t139;
t185 = pkin(2) * t200;
t232 = pkin(3) * t129;
t82 = -t185 - t232;
t234 = -t136 * t82 + t139 * t159;
t140 = -pkin(8) - pkin(7);
t100 = g(1) * t106;
t115 = sin(t127);
t231 = g(1) * t115;
t229 = g(3) * t138;
t208 = t134 * t135;
t183 = t129 * t208;
t205 = t137 * t138;
t57 = -t129 * t205 + t183;
t227 = t59 * t57;
t111 = pkin(2) * t136 + pkin(7);
t226 = -pkin(8) - t111;
t94 = t140 * t135;
t123 = t138 * pkin(8);
t95 = pkin(7) * t138 + t123;
t47 = -t134 * t95 + t137 * t94;
t69 = -t205 + t208;
t181 = qJD(4) * t140;
t78 = t135 * t181;
t79 = t138 * t181;
t225 = qJD(5) * t47 + t134 * t79 + t137 * t78 + t185 * t69;
t48 = t134 * t94 + t137 * t95;
t224 = -qJD(5) * t48 - t134 * t78 + t137 * t79 + t185 * t70;
t207 = t135 * t126;
t165 = t134 * t207 - t137 * t204;
t43 = t128 * t70;
t21 = t129 * t43 + t165;
t160 = t128 * t208;
t42 = t160 + t237;
t223 = -t21 * t70 + t42 * t57;
t222 = t100 * t138 + t197 * t82;
t221 = pkin(3) * t107 + pkin(7) * t106;
t219 = t134 * t46;
t216 = t137 * t46;
t133 = qJ(4) + qJ(5);
t121 = sin(t133);
t213 = t106 * t121;
t122 = cos(t133);
t212 = t106 * t122;
t211 = t107 * t121;
t210 = t107 * t122;
t209 = t129 * t135;
t202 = t130 - t131;
t190 = t135 * t238 + t82 * t196;
t187 = pkin(2) * t198;
t186 = pkin(4) * t197;
t124 = t129 ^ 2;
t182 = t135 * t124 * t138;
t112 = pkin(4) * t138 + pkin(3);
t180 = t129 * t199;
t173 = qJD(4) * t226;
t172 = -t126 * t206 + t129 * t237 - t134 * t204;
t170 = -t106 * t140 + t107 * t112;
t169 = t201 * t126;
t168 = t129 * t189;
t167 = t138 * t179;
t166 = g(1) * (-pkin(3) * t106 + pkin(7) * t107);
t116 = cos(t127);
t163 = -g(2) * t116 + t231;
t20 = t129 * t160 + t172;
t162 = -t20 * t69 + t43 * t59;
t18 = t137 * t41 - t219;
t19 = t134 * t41 + t216;
t5 = -qJD(5) * t19 - t134 * t17 + t137 * t16;
t161 = t18 * t42 - t19 * t43 - t4 * t69 - t5 * t70 - t236;
t125 = qJDD(4) + qJDD(5);
t28 = t125 * t70 - t128 * t42;
t67 = t226 * t135;
t68 = t111 * t138 + t123;
t38 = -t134 * t68 + t137 * t67;
t39 = t134 * t67 + t137 * t68;
t158 = -t106 * t112 - t107 * t140;
t157 = t100 + t203 - t230;
t37 = pkin(4) * t154 + t61;
t60 = -t112 * t129 - t185;
t156 = -g(1) * t213 + g(2) * t211 + t37 * t70 - t42 * t60;
t155 = g(1) * t212 - g(2) * t210 + t37 * t69 + t43 * t60;
t153 = t186 - t189;
t152 = -t129 * t82 + t236;
t141 = qJD(4) ^ 2;
t151 = pkin(7) * t141 - t168 - t228;
t113 = -pkin(3) - t233;
t149 = pkin(2) * t180 + t111 * t141 + t113 * t126;
t148 = -pkin(7) * qJDD(4) + (t185 - t232) * qJD(4);
t147 = -qJDD(4) * t111 + (t113 * t129 - t187) * qJD(4);
t145 = -t236 + t146;
t144 = g(1) * t210 + g(2) * t212 + g(3) * t121 + t57 * t60 - t4;
t143 = g(1) * t211 + g(2) * t213 - g(3) * t122 - t59 * t60 + t5;
t132 = qJDD(1) - g(3);
t105 = pkin(2) * t116;
t89 = -t112 - t233;
t86 = qJDD(4) * t138 - t135 * t141;
t85 = qJDD(4) * t135 + t138 * t141;
t80 = t186 + t188;
t64 = t126 * t131 - 0.2e1 * t167;
t63 = t126 * t130 + 0.2e1 * t167;
t50 = -t135 * t187 + t138 * t173;
t49 = t135 * t173 + t138 * t187;
t44 = -0.2e1 * qJD(4) * t129 * t202 + 0.2e1 * t135 * t204;
t29 = -t125 * t69 - t128 * t43;
t24 = -t57 ^ 2 + t59 ^ 2;
t23 = t137 * t45 - t219;
t22 = -t134 * t45 - t216;
t12 = -t172 + (-t183 + t57) * t128;
t11 = -qJD(5) * t39 - t134 * t49 + t137 * t50;
t10 = qJD(5) * t38 + t134 * t50 + t137 * t49;
t7 = t21 * t69 + t43 * t57;
t6 = -t20 * t70 - t42 * t59;
t1 = -t162 + t223;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t132, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, 0, 0, 0, 0, 0, 0, t86, -t85, 0, -qJD(4) * t159 + t31 * t135 + t32 * t138 - g(3), 0, 0, 0, 0, 0, 0, t29, -t28, t162 + t223, -t18 * t43 - t19 * t42 + t4 * t70 - t5 * t69 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t163, g(1) * t116 + g(2) * t115, 0, 0, 0, 0, 0, 0, 0, t126, (t126 * t139 - t180) * pkin(2) + t157, ((-qJDD(2) - t126) * t136 + (-qJD(2) - t129) * t198) * pkin(2) + t236, 0, (t163 + (t136 ^ 2 + t139 ^ 2) * qJDD(2) * pkin(2)) * pkin(2), t63, t44, t85, t64, t86, 0, t147 * t135 + (-t149 - t238) * t138 + t222, t147 * t138 + (t149 - t100) * t135 + t190, t111 * t169 + t187 * t235 + t145, t61 * t113 - t166 - g(2) * (t105 + t221) + (-qJD(3) * t234 + t231) * pkin(2) + t146 * t111, t6, t1, t28, t7, t29, 0, t11 * t128 + t125 * t38 + t21 * t89 + t57 * t80 + t155, -t10 * t128 - t125 * t39 - t20 * t89 + t59 * t80 + t156, -t10 * t57 - t11 * t59 + t20 * t38 - t21 * t39 + t161, t4 * t39 + t19 * t10 + t5 * t38 + t18 * t11 + t37 * t89 + t60 * t80 - g(1) * (-pkin(2) * t115 + t158) - g(2) * (t105 + t170); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, t157 + t168, (-t191 + (-qJD(3) + t129) * t200) * pkin(2) + t236, 0, 0, t63, t44, t85, t64, t86, 0, t148 * t135 + (-t151 - t238) * t138 + t222, t148 * t138 + (t151 - t100) * t135 + t190, pkin(7) * t169 - t185 * t235 + t145, -t61 * pkin(3) + pkin(7) * t146 - g(2) * t221 + t220 * t234 - t166, t6, t1, t28, t7, t29, 0, -t112 * t21 + t125 * t47 + t128 * t224 + t153 * t57 + t155, t112 * t20 - t125 * t48 - t128 * t225 + t153 * t59 + t156, t20 * t47 - t21 * t48 - t224 * t59 - t225 * t57 + t161, -g(1) * t158 - g(2) * t170 - t37 * t112 + t153 * t60 + t18 * t224 + t19 * t225 + t4 * t48 + t5 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t202 * t124, t207, t182, t204, qJDD(4), -t229 + t118 + (t56 - t215) * qJD(4) + (t152 - t239) * t135, g(3) * t135 + (t55 + t218) * qJD(4) + t152 * t138 - t184, 0, 0, t227, t24, t12, -t227, -t165, t125, -t128 * t22 + (t125 * t137 - t128 * t195 - t209 * t57) * pkin(4) + t143, t128 * t23 + (-t125 * t134 - t128 * t194 - t209 * t59) * pkin(4) + t144, (t19 + t22) * t59 + (-t18 + t23) * t57 + (-t134 * t21 + t137 * t20 + (t134 * t59 - t137 * t57) * qJD(5)) * pkin(4), -t18 * t22 - t19 * t23 + (-t229 + t134 * t4 + t137 * t5 + (-t134 * t18 + t137 * t19) * qJD(5) + (-t129 * t60 + t236) * t135) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t24, t12, -t227, -t165, t125, t128 * t19 + t143, t128 * t18 + t144, 0, 0;];
tau_reg = t2;
