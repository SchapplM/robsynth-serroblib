% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRP5
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:25
% EndTime: 2019-12-05 16:49:29
% DurationCPUTime: 2.30s
% Computational Cost: add. (2457->321), mult. (5599->396), div. (0->0), fcn. (3965->10), ass. (0->184)
t128 = sin(qJ(4));
t129 = sin(qJ(3));
t131 = cos(qJ(3));
t234 = cos(qJ(4));
t83 = t128 * t131 + t234 * t129;
t73 = t83 * qJD(2);
t215 = t73 * qJ(5);
t130 = sin(qJ(2));
t197 = t130 * qJD(1);
t102 = qJD(2) * pkin(6) + t197;
t172 = pkin(7) * qJD(2) + t102;
t66 = t172 * t131;
t53 = t128 * t66;
t65 = t172 * t129;
t56 = qJD(3) * pkin(3) - t65;
t35 = t234 * t56 - t53;
t17 = t35 - t215;
t174 = qJDD(2) * t234;
t192 = t129 * qJDD(2);
t122 = qJD(3) + qJD(4);
t247 = t122 * t83;
t31 = qJD(2) * t247 + t128 * t192 - t131 * t174;
t132 = cos(qJ(2));
t126 = sin(pkin(8));
t127 = cos(pkin(8));
t164 = g(1) * t127 + g(2) * t126;
t156 = t164 * t132;
t230 = g(3) * t130;
t147 = t156 + t230;
t157 = t164 * t130;
t205 = qJDD(1) - g(3);
t246 = t205 * t132 + t157;
t196 = t132 * qJD(1);
t219 = qJD(2) * pkin(2);
t103 = -t196 - t219;
t123 = t129 ^ 2;
t124 = t131 ^ 2;
t203 = t123 + t124;
t171 = t203 * t132;
t245 = t102 * t171 + t103 * t130;
t182 = t234 * t131;
t207 = t128 * t129;
t153 = t182 - t207;
t150 = t153 * t132;
t179 = t234 * qJD(4);
t198 = qJD(4) * t128;
t237 = pkin(7) + pkin(6);
t183 = qJD(3) * t237;
t88 = t129 * t183;
t89 = t131 * t183;
t94 = t237 * t129;
t95 = t237 * t131;
t223 = -qJD(1) * t150 - t128 * t89 - t94 * t179 - t95 * t198 - t234 * t88;
t49 = -t128 * t94 + t234 * t95;
t222 = -t49 * qJD(4) + t128 * t88 + t83 * t196 - t234 * t89;
t120 = qJDD(3) + qJDD(4);
t167 = pkin(3) * t179;
t233 = pkin(3) * t128;
t244 = -t120 * t233 - t122 * t167;
t68 = t153 * t130;
t117 = t120 * pkin(4);
t162 = t122 * t207;
t166 = qJD(2) * t182;
t191 = t131 * qJDD(2);
t169 = -t122 * t166 - t128 * t191 - t129 * t174;
t30 = qJD(2) * t162 + t169;
t218 = t30 * qJ(5);
t241 = t117 + t218;
t125 = qJ(3) + qJ(4);
t118 = sin(t125);
t208 = t127 * t132;
t119 = cos(t125);
t209 = t127 * t119;
t210 = t126 * t132;
t211 = t126 * t119;
t240 = t118 * t230 - g(2) * (-t118 * t210 - t209) - g(1) * (-t118 * t208 + t211);
t134 = qJD(3) ^ 2;
t229 = g(3) * t132;
t195 = qJD(1) * qJD(2);
t112 = t130 * t195;
t190 = t132 * qJDD(1);
t213 = qJDD(2) * pkin(2);
t78 = t112 - t190 - t213;
t180 = -t78 - t229;
t239 = -pkin(6) * t134 + t130 * (t164 + t195) + t180 + t213;
t238 = t73 ^ 2;
t228 = t131 * pkin(3);
t201 = qJD(2) * t129;
t181 = t128 * t201;
t71 = -t166 + t181;
t227 = t73 * t71;
t226 = -qJ(5) * t247 + qJD(5) * t153 + t223;
t44 = -qJD(3) * t182 - t131 * t179 + t162;
t225 = t44 * qJ(5) - t83 * qJD(5) + t222;
t16 = t122 * pkin(4) + t17;
t224 = -t17 + t16;
t221 = -t71 * t167 - t31 * t233;
t39 = -t234 * t65 - t53;
t220 = qJ(5) * t31;
t217 = t71 * qJ(5);
t216 = t71 * t122;
t214 = t73 * t122;
t176 = t71 * pkin(4) + qJD(5);
t116 = pkin(2) + t228;
t80 = -t116 * qJD(2) - t196;
t43 = t176 + t80;
t206 = qJD(5) + t43;
t204 = t123 - t124;
t135 = qJD(2) ^ 2;
t202 = t134 + t135;
t200 = qJD(2) * t130;
t199 = qJD(3) * t129;
t194 = qJD(2) * qJD(3);
t193 = qJDD(3) * t129;
t189 = t132 * qJDD(2);
t188 = t234 * pkin(3);
t187 = pkin(3) * t198;
t186 = pkin(3) * t199;
t185 = pkin(3) * t201;
t55 = t234 * t66;
t184 = t129 * t135 * t131;
t178 = t129 * t194;
t177 = t131 * t194;
t93 = pkin(4) * t119 + t228;
t79 = qJDD(2) * pkin(6) + t130 * qJDD(1) + t132 * t195;
t34 = -t131 * t102 * qJD(3) + qJDD(3) * pkin(3) - t129 * t79 + (-t177 - t192) * pkin(7);
t37 = -t102 * t199 + t131 * t79 + (-t178 + t191) * pkin(7);
t175 = -t128 * t37 + t234 * t34;
t38 = t128 * t65 - t55;
t48 = -t128 * t95 - t234 * t94;
t173 = t203 * t79;
t7 = t128 * t34 + t56 * t179 - t66 * t198 + t234 * t37;
t42 = pkin(4) * t247 + t186;
t170 = t42 - t197;
t168 = t203 * qJDD(2);
t165 = t129 * t177;
t163 = g(1) * t126 - g(2) * t127;
t159 = -t119 * t229 + (g(1) * t209 + g(2) * t211) * t130;
t36 = t128 * t56 + t55;
t155 = t163 * t131;
t154 = t186 - t197;
t152 = pkin(3) * t178 - t116 * qJDD(2) + t112;
t148 = -g(1) * (-t126 * t118 - t119 * t208) - g(2) * (t127 * t118 - t119 * t210) + t119 * t230 - t7;
t145 = -t122 * t181 - t169;
t144 = t31 * pkin(4) + qJDD(5) + t152;
t143 = -pkin(6) * qJDD(3) + (t103 + t196 - t219) * qJD(3);
t142 = t80 * t71 + t148;
t8 = -t36 * qJD(4) + t175;
t141 = t118 * t229 + (-qJD(1) * t73 - t164 * t118) * t130;
t140 = -t103 * qJD(2) + t147 - t79;
t139 = t206 * t71 + t148 + t220;
t138 = t8 + t240;
t136 = -t80 * t73 + t138;
t121 = -qJ(5) - t237;
t115 = t188 + pkin(4);
t92 = -t129 * pkin(3) - pkin(4) * t118;
t87 = pkin(2) + t93;
t70 = t71 ^ 2;
t67 = t83 * t130;
t60 = -pkin(4) * t153 - t116;
t51 = t73 * pkin(4) + t185;
t46 = t152 - t190;
t41 = qJ(5) * t153 + t49;
t40 = -t83 * qJ(5) + t48;
t32 = -t70 + t238;
t27 = t120 * t153 - t122 * t247;
t26 = t120 * t83 - t122 * t44;
t22 = -t122 * t68 - t132 * t73;
t21 = qJD(2) * t150 - t130 * t247;
t20 = -t215 + t39;
t19 = t38 + t217;
t18 = t36 - t217;
t15 = t214 - t31;
t14 = t145 + t216;
t13 = t144 - t190;
t10 = -t153 * t31 + t247 * t71;
t9 = -t30 * t83 - t44 * t73;
t6 = -t67 * t120 + t22 * t122 - t132 * t31 + t71 * t200;
t5 = -t68 * t120 - t21 * t122 + t132 * t30 + t73 * t200;
t4 = -qJD(5) * t71 - t220 + t7;
t3 = -t153 * t30 - t247 * t73 - t31 * t83 + t44 * t71;
t2 = -t73 * qJD(5) + t241 + t8;
t1 = -t21 * t71 - t22 * t73 - t30 * t67 - t31 * t68;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t205, 0, 0, 0, 0, 0, 0, -t135 * t130 + t189, -qJDD(2) * t130 - t135 * t132, 0, -g(3) + (t130 ^ 2 + t132 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, (-0.2e1 * t178 + t191) * t132 + (-t202 * t131 - t193) * t130, (-qJDD(3) * t130 - 0.2e1 * t132 * t194) * t131 + (t202 * t130 - t189) * t129, t130 * t168 + t135 * t171, t245 * qJD(2) + t130 * t173 - t78 * t132 - g(3), 0, 0, 0, 0, 0, 0, t6, t5, t1, -t46 * t132 + t80 * t200 + t36 * t21 + t35 * t22 - t8 * t67 + t7 * t68 - g(3), 0, 0, 0, 0, 0, 0, t6, t5, t1, -t13 * t132 + t16 * t22 + t18 * t21 - t2 * t67 + t43 * t200 + t4 * t68 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t246, -t205 * t130 + t156, 0, 0, t123 * qJDD(2) + 0.2e1 * t165, 0.2e1 * t129 * t191 - 0.2e1 * t204 * t194, t134 * t131 + t193, t124 * qJDD(2) - 0.2e1 * t165, qJDD(3) * t131 - t134 * t129, 0, t143 * t129 + t239 * t131, -t239 * t129 + t143 * t131, -t230 + t173 + pkin(6) * t168 + (-t203 * t195 - t164) * t132, (t157 + t180) * pkin(2) + (t173 - t147) * pkin(6) - t245 * qJD(1), t9, t3, t26, t10, t27, 0, -t116 * t31 + t48 * t120 + t222 * t122 - t153 * t46 + t154 * t71 + t247 * t80 + t159, t116 * t30 - t49 * t120 - t223 * t122 + t73 * t186 - t80 * t44 + t46 * t83 + t141, t153 * t7 - t222 * t73 - t223 * t71 - t247 * t36 + t48 * t30 - t49 * t31 + t35 * t44 - t8 * t83 - t147, t7 * t49 + t8 * t48 - t46 * t116 - g(3) * (t132 * t116 + t130 * t237) + t154 * t80 + t223 * t36 + t222 * t35 + t164 * (t116 * t130 - t132 * t237), t9, t3, t26, t10, t27, 0, t40 * t120 + t225 * t122 - t13 * t153 + t170 * t71 + t247 * t43 + t60 * t31 + t159, -t41 * t120 - t226 * t122 + t13 * t83 - t60 * t30 + t42 * t73 - t43 * t44 + t141, t153 * t4 + t16 * t44 - t18 * t247 - t2 * t83 - t225 * t73 - t226 * t71 + t40 * t30 - t41 * t31 - t147, t4 * t41 + t2 * t40 + t13 * t60 - g(3) * (-t130 * t121 + t132 * t87) + t170 * t43 + t226 * t18 + t225 * t16 + t164 * (t121 * t132 + t130 * t87); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, t204 * t135, t192, t184, t191, qJDD(3), t129 * t140 - t155, t163 * t129 + t140 * t131, 0, 0, t227, t32, t14, -t227, t15, t120, -t38 * t122 + (t234 * t120 - t122 * t198 - t71 * t201) * pkin(3) + t136, t39 * t122 - t185 * t73 + t142 + t244, t30 * t188 + (-t35 + t39) * t71 + (t36 + t38 + t187) * t73 + t221, -t35 * t38 - t36 * t39 + (t8 * t234 + t7 * t128 - t155 + (-t80 * qJD(2) + t147) * t129 + (-t35 * t128 + t36 * t234) * qJD(4)) * pkin(3), t227, t32, t14, -t227, t15, t120, t115 * t120 - t19 * t122 - t51 * t71 - t206 * t73 + (-t55 + (-pkin(3) * t122 - t56) * t128) * qJD(4) + t175 + t240 + t241, t122 * t20 - t51 * t73 + t139 + t244, t115 * t30 + (-t16 + t20) * t71 + (t18 + t19 + t187) * t73 + t221, t2 * t115 - t18 * t20 - t16 * t19 - t43 * t51 - g(1) * (t126 * t93 + t92 * t208) - g(2) * (-t127 * t93 + t92 * t210) - t92 * t230 + (t4 * t128 + (-t128 * t16 + t234 * t18) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t32, t14, -t227, t15, t120, t36 * t122 + t136, t122 * t35 + t142, 0, 0, t227, t32, t14, -t227, t15, t120, t218 + t18 * t122 + 0.2e1 * t117 + (-t176 - t43) * t73 + t138, -pkin(4) * t238 + t122 * t17 + t139, pkin(4) * t30 - t224 * t71, t224 * t18 + (-t43 * t73 + t2 + t240) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 + t214, t145 - t216, -t70 - t238, t16 * t73 + t18 * t71 + t144 - t246;];
tau_reg = t11;
