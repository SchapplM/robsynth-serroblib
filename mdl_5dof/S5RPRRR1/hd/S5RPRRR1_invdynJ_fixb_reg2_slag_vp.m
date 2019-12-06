% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:09:58
% EndTime: 2019-12-05 18:10:09
% DurationCPUTime: 3.57s
% Computational Cost: add. (1841->381), mult. (5015->529), div. (0->0), fcn. (3776->8), ass. (0->208)
t80 = cos(qJ(4));
t196 = qJD(4) * t80;
t81 = cos(qJ(3));
t214 = t80 * t81;
t75 = sin(qJ(5));
t77 = sin(qJ(3));
t222 = t75 * t77;
t79 = cos(qJ(5));
t110 = t79 * t214 + t222;
t29 = t110 * qJD(1);
t267 = t79 * t196 - t29;
t70 = t81 * qJD(1);
t250 = t70 - qJD(4);
t249 = qJD(4) * t250;
t125 = t80 * t249;
t67 = t81 * qJDD(1);
t181 = qJDD(4) - t67;
t76 = sin(qJ(4));
t89 = -t76 * t181 + t125;
t126 = t76 * t249;
t266 = t80 * t181 + t126;
t177 = qJD(1) * qJD(4);
t158 = t80 * t177;
t142 = qJD(4) + t70;
t64 = t77 * qJDD(1);
t90 = t142 * qJD(3) + t64;
t16 = -t80 * qJDD(3) + t77 * t158 + t90 * t76;
t15 = -t80 * t90 + (t77 * t177 - qJDD(3)) * t76;
t185 = qJ(2) * qJD(1);
t163 = t81 * t185;
t187 = t80 * qJD(2);
t44 = t76 * t163 - t187;
t116 = t250 * t44;
t189 = t76 * qJD(2);
t45 = t80 * t163 + t189;
t265 = t45 * t250;
t255 = t81 * t250;
t264 = t80 * t255;
t193 = qJD(5) * t75;
t263 = t76 * t193 - t267;
t73 = t77 ^ 2;
t262 = t73 * qJD(1) + t255;
t78 = sin(qJ(1));
t82 = cos(qJ(1));
t260 = g(1) * t78 - g(2) * t82;
t162 = qJDD(2) - t260;
t200 = qJD(1) * t77;
t68 = t80 * qJD(3);
t46 = t76 * t200 - t68;
t42 = qJD(5) + t46;
t261 = t263 * t42;
t153 = t79 * t250;
t166 = t80 * t200;
t188 = t76 * qJD(3);
t48 = t166 + t188;
t21 = t48 * t75 + t153;
t258 = t250 * t21;
t178 = qJD(1) * qJD(3);
t159 = t77 * t178;
t254 = (t159 - t67) * qJ(2);
t176 = qJD(2) * qJD(1);
t180 = qJ(2) * qJDD(1);
t253 = t176 + t180;
t164 = t77 * t185;
t24 = t79 * t164 - t75 * t45;
t252 = -t44 * t21 - t24 * t42;
t243 = g(2) * t78;
t245 = g(1) * t82;
t251 = -t245 - t243;
t114 = t142 * qJD(2);
t156 = t81 * t178;
t136 = qJ(2) * t156;
t157 = t77 * t176;
t218 = t77 * t79;
t25 = t75 * t164 + t79 * t45;
t135 = qJD(4) * t163;
t8 = (qJDD(2) - t135) * t76 + (t114 - t254) * t80;
t2 = -qJD(5) * t25 + t180 * t218 - t75 * t8 + (t136 + t157) * t79;
t23 = -t250 * t75 + t79 * t48;
t248 = t44 * t23 - t25 * t42 - t2;
t223 = t75 * t76;
t105 = t159 + t181;
t4 = qJD(5) * t153 - t75 * t105 + t79 * t15 + t48 * t193;
t212 = t81 * t75;
t174 = t80 * t212;
t28 = qJD(1) * t174 - t79 * t200;
t128 = t75 * t196 - t28;
t192 = qJD(5) * t79;
t98 = t76 * t192 + t128;
t247 = t4 * t223 - t23 * t98;
t140 = -qJD(5) + t68;
t195 = qJD(4) * t81;
t171 = t76 * t195;
t246 = t140 * t77 + t171;
t5 = qJD(5) * t23 - t79 * t105 - t75 * t15;
t210 = t82 * t76;
t39 = t78 * t214 - t210;
t244 = g(2) * t39;
t242 = g(3) * t77;
t241 = g(3) * t81;
t240 = t5 * t80;
t14 = qJDD(5) + t16;
t239 = t14 * t75;
t238 = t14 * t79;
t237 = t15 * t76;
t236 = t16 * t80;
t235 = t21 * t42;
t234 = t23 * t21;
t233 = t23 * t42;
t228 = t46 * t76;
t227 = t46 * t80;
t226 = t48 * t46;
t225 = t48 * t76;
t224 = t48 * t80;
t221 = t76 * t77;
t220 = t76 * t79;
t219 = t76 * t81;
t217 = t77 * t80;
t216 = t77 * t82;
t215 = t79 * t81;
t213 = t80 * t82;
t84 = qJD(1) ^ 2;
t211 = t81 * t84;
t209 = -t76 * t16 - t46 * t196;
t138 = t76 * t164;
t208 = -qJD(3) * t138 - t80 * qJDD(2);
t165 = 0.2e1 * t176;
t182 = t73 * qJDD(1);
t85 = qJ(2) ^ 2;
t207 = t73 * qJ(2) * t165 + t85 * t182;
t74 = t81 ^ 2;
t205 = t73 - t74;
t204 = t73 + t74;
t203 = qJ(2) * t77;
t202 = qJ(2) * t81;
t201 = qJ(2) * t84;
t199 = qJD(3) * t77;
t198 = qJD(3) * t81;
t197 = qJD(4) * t76;
t194 = qJD(5) * t45;
t191 = t46 * qJD(3);
t184 = qJ(2) * qJD(3);
t183 = qJDD(1) * t85;
t179 = qJ(2) * qJDD(3);
t173 = t77 * t211;
t172 = t23 * t70;
t170 = t80 * t195;
t169 = t77 * t197;
t168 = t46 * t200;
t161 = t23 * t197 + t4 * t80;
t160 = t204 * t84;
t155 = t81 * t180;
t150 = t42 * t79;
t149 = -t48 + t188;
t148 = qJD(1) * t250;
t147 = qJD(3) * t250;
t144 = 0.2e1 * t159;
t141 = -qJD(5) * t80 + qJD(3);
t137 = t77 * t156;
t38 = t78 * t219 + t213;
t40 = t81 * t210 - t78 * t80;
t133 = g(1) * t40 + g(2) * t38;
t132 = -g(1) * t38 + g(2) * t40;
t41 = t81 * t213 + t76 * t78;
t131 = -g(1) * t41 - t244;
t124 = -t24 * t79 - t25 * t75;
t123 = t24 * t75 - t25 * t79;
t122 = t44 * t80 - t45 * t76;
t121 = t44 * t76 + t45 * t80;
t120 = -t225 + t227;
t119 = t225 + t227;
t117 = t76 * t250;
t115 = t141 * t81;
t113 = qJD(5) * t164 + t8;
t111 = t78 * t218 - t39 * t75;
t37 = t79 * t217 - t212;
t109 = -t174 + t218;
t36 = t75 * t217 + t215;
t108 = t128 * t42;
t107 = -t15 * t80 - t48 * t197;
t106 = -t42 * t192 - t239;
t104 = qJD(3) * (t142 - t70);
t103 = qJ(2) * (qJD(1) * t70 - t211);
t101 = t80 * t105 - t148 * t219 + t126;
t100 = t81 * t188 + t77 * t196;
t99 = t251 + t165;
t96 = t263 * t21 - t5 * t220;
t95 = t243 - t253;
t83 = qJD(3) ^ 2;
t94 = -qJ(2) * t83 - t162;
t9 = t80 * t135 + (t114 + t155) * t76 + t208;
t93 = g(3) * t221 + t133 - t9;
t92 = t242 - t114;
t88 = -t99 - t180;
t87 = t95 + t245;
t86 = t122 * qJD(4) + t9 * t76 + t8 * t80;
t31 = t110 * qJ(2);
t30 = t109 * qJ(2);
t27 = t37 * t185;
t26 = t36 * t185;
t20 = t75 * t216 + t41 * t79;
t19 = t79 * t216 - t41 * t75;
t11 = t140 * t215 + (t141 * t75 - t79 * t197) * t77;
t10 = -t75 * t169 - t81 * t193 - t79 * t199 + (t77 * t192 + t75 * t198) * t80;
t7 = t109 * qJD(2) + (t79 * t115 + t246 * t75) * qJ(2);
t6 = t110 * qJD(2) + (t75 * t115 - t246 * t79) * qJ(2);
t1 = t113 * t79 + (t157 - t194 + (t156 + t64) * qJ(2)) * t75;
t3 = [0, 0, 0, 0, 0, qJDD(1), t260, -t251, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -t162, 0, t99 + 0.2e1 * t180, qJ(2) * t99 + t183, 0.2e1 * t137 + t182, -0.2e1 * t205 * t178 + 0.2e1 * t77 * t67, qJDD(3) * t77 + t81 * t83, qJDD(1) * t74 - 0.2e1 * t137, qJDD(3) * t81 - t77 * t83, 0, -t77 * t179 + t81 * t94, -t81 * t179 - t77 * t94, 0.2e1 * t253 * t204 + t251, t74 * t183 + (t74 * t165 + t251) * qJ(2) + t207, -t15 * t217 + (t81 * t68 - t169) * t48, -t119 * t198 + (t237 - t236 + (-t224 + t228) * qJD(4)) * t77, (-t147 * t80 + t15) * t81 + ((t48 + t166) * qJD(3) + t266) * t77, t100 * t46 + t16 * t221, t16 * t81 + (t125 - t191) * t77 + (-t105 * t77 + t147 * t81) * t76, -t105 * t81 - t147 * t77, g(1) * t39 - g(2) * t41 + t9 * t81 + t262 * t189 + ((t76 * qJDD(1) + t158) * t73 + (t89 + t191) * t81) * qJ(2) + (qJD(2) * t46 - t44 * qJD(3) + (t104 * t76 + t16) * qJ(2)) * t77, t8 * t81 + t262 * t187 + ((t80 * qJDD(1) - t177 * t76) * t73 + (t48 * qJD(3) - t266) * t81) * qJ(2) + (qJD(2) * t48 - t45 * qJD(3) + (t104 * t80 - t15) * qJ(2)) * t77 + t132, (t122 * qJD(3) - t120 * qJD(2) + (-t237 - t236 + (t224 + t228) * qJD(4)) * qJ(2)) * t81 + ((-qJD(4) * t45 + t184 * t46 + t9) * t80 + (-qJD(4) * t44 - t184 * t48 - t8) * t76 + t260) * t77, (qJD(2) * t121 + t144 * t85) * t81 + (-t121 * t199 + t81 * t86 + t251) * qJ(2) + t207, t11 * t23 - t37 * t4, -t10 * t23 - t11 * t21 + t36 * t4 - t37 * t5, t100 * t23 + t11 * t42 + t14 * t37 - t4 * t221, t10 * t21 + t36 * t5, -t10 * t42 - t100 * t21 - t14 * t36 - t5 * t221, t100 * t42 + t14 * t221, t7 * t42 + t30 * t14 + t9 * t36 + t44 * t10 - g(1) * (-t78 * t222 - t39 * t79) - g(2) * t20 + (t202 * t21 + t24 * t77) * t196 + ((-t184 * t21 + t2) * t77 + (qJ(2) * t5 + qJD(2) * t21 + qJD(3) * t24) * t81) * t76, -t6 * t42 - t31 * t14 + t9 * t37 + t44 * t11 + g(1) * t111 - g(2) * t19 + (t202 * t23 - t25 * t77) * t196 + ((-t184 * t23 - t1) * t77 + (-qJ(2) * t4 + qJD(2) * t23 - qJD(3) * t25) * t81) * t76, -t1 * t36 - t10 * t25 - t11 * t24 - t2 * t37 - t21 * t6 - t23 * t7 + t30 * t4 - t31 * t5 - t132, t44 * t81 * t189 + t1 * t31 + t2 * t30 + t24 * t7 + t25 * t6 + (t9 * t219 + (-t188 * t77 + t170) * t44 + t251) * qJ(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t84, t162 - t201, 0, 0, 0, 0, 0, 0, -t67 + t144, t64 + 0.2e1 * t156, -t160, -qJ(2) * t160 + t162, 0, 0, 0, 0, 0, 0, t101 - t168, (-t264 + (-t48 - t188) * t77) * qJD(1) + t89, t120 * t70 - t107 + t209, -t73 * t201 + (-t9 - t265) * t80 + (t8 - t116) * t76 - t260, 0, 0, 0, 0, 0, 0, -t240 - t108 + (t106 - t258) * t76, (-t172 - t238) * t76 + t261 + t161, -t247 + t96, t24 * t28 - t25 * t29 + (-qJD(4) * t123 - t9) * t80 + (qJD(5) * t124 + t1 * t79 - t2 * t75 - t116) * t76 - t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, t205 * t84, t64, t173, t67, qJDD(3), t77 * t88 - t241, t81 * t88 + t242, 0, 0, -t224 * t250 - t237, t119 * t70 + t107 + t209, (t149 * t77 + t264) * qJD(1) - t89, -t117 * t46 - t236, t101 + t168, t77 * t148, (-g(3) * t80 + (-t46 - t68) * t185) * t81 + (t44 * qJD(1) + t103 * t76 + t80 * t87) * t77, (g(3) * t76 + t149 * t185) * t81 + (t45 * qJD(1) + t103 * t80 - t76 * t87) * t77, -t242 + t251 * t81 + (-t120 * t203 - t122 * t81) * qJD(1) + t86, (t121 * t185 - t211 * t85) * t77, -t4 * t220 - t23 * t263, t247 + t96, (-t172 + t238) * t76 - t261 + t161, t21 * t98 + t5 * t223, t240 - t108 + (t106 + t258) * t76, -t117 * t42 - t14 * t80, -t2 * t80 - t26 * t42 - g(3) * t110 + t128 * t44 + (t44 * t192 + t24 * qJD(4) + t9 * t75 + (t203 * t21 - t24 * t81) * qJD(1)) * t76 - t251 * t37, t1 * t80 - t27 * t42 - g(3) * t109 + t267 * t44 + (-t44 * t193 - t25 * qJD(4) + t9 * t79 + (t203 * t23 + t25 * t81) * qJD(1)) * t76 + t251 * t36, -t21 * t27 + t23 * t26 + t24 * t29 + t25 * t28 + t124 * t196 + (qJD(5) * t123 - t1 * t75 - t2 * t79 - t251 * t77 - t241) * t76, t138 * t44 - t24 * t26 + t25 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, -t46 ^ 2 + t48 ^ 2, -t250 * t46 - t15, -t226, -t250 * t48 - t16, t105, -t265 + (-t48 * t77 - t170) * t185 + (t92 - t155) * t76 + t133 - t208, -t76 * qJDD(2) + t116 + (t46 * t77 + t171) * t185 + (t92 + t254) * t80 - t131, 0, 0, t150 * t23 - t4 * t75, (-t4 - t235) * t79 + (-t5 - t233) * t75, t150 * t42 - t23 * t48 + t239, t75 * t235 - t5 * t79, -t42 ^ 2 * t75 + t21 * t48 + t238, -t42 * t48, -t21 * t45 - t24 * t48 + t79 * t93, -t23 * t45 + t25 * t48 - t75 * t93, -g(3) * t217 + (t1 + t252) * t79 + t248 * t75 + t131, (-t123 - t45) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, -t21 ^ 2 + t23 ^ 2, -t4 + t235, -t234, t233 - t5, t14, -g(1) * t19 - g(2) * t111 + g(3) * t36 - t248, g(1) * t20 + g(3) * t37 + (-t113 + t244) * t79 + (t77 * t95 - t136 + t194) * t75 - t252, 0, 0;];
tau_reg = t3;
