% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:30
% EndTime: 2019-12-05 16:23:37
% DurationCPUTime: 2.93s
% Computational Cost: add. (3279->352), mult. (7587->468), div. (0->0), fcn. (5688->14), ass. (0->181)
t153 = sin(qJ(3));
t156 = cos(qJ(3));
t236 = qJ(4) + pkin(6);
t194 = qJD(3) * t236;
t100 = -t153 * qJD(4) - t156 * t194;
t147 = sin(pkin(9));
t149 = cos(pkin(9));
t113 = t147 * t156 + t149 * t153;
t157 = cos(qJ(2));
t210 = t157 * qJD(1);
t99 = t156 * qJD(4) - t153 * t194;
t235 = t149 * t100 + t113 * t210 - t147 * t99;
t222 = t149 * t156;
t112 = t147 * t153 - t222;
t174 = t112 * t157;
t234 = qJD(1) * t174 + t147 * t100 + t149 * t99;
t152 = sin(qJ(5));
t155 = cos(qJ(5));
t200 = qJD(2) * t222;
t216 = qJD(2) * t153;
t103 = t147 * t216 - t200;
t106 = t113 * qJD(2);
t183 = t152 * t103 - t155 * t106;
t105 = t113 * qJD(3);
t205 = t156 * qJDD(2);
t206 = t153 * qJDD(2);
t184 = t147 * t206 - t149 * t205;
t57 = qJD(2) * t105 + t184;
t208 = qJD(2) * qJD(3);
t198 = t153 * t208;
t170 = qJDD(2) * t113 - t147 * t198;
t197 = t156 * t208;
t58 = t149 * t197 + t170;
t165 = qJD(5) * t183 - t152 * t58 - t155 * t57;
t143 = qJD(3) + qJD(5);
t228 = t183 * t143;
t263 = t165 - t228;
t213 = qJD(5) * t155;
t214 = qJD(5) * t152;
t175 = -t103 * t213 - t106 * t214 - t152 * t57 + t155 * t58;
t52 = -t155 * t103 - t152 * t106;
t227 = t52 * t143;
t262 = t175 - t227;
t238 = t183 ^ 2;
t239 = t52 ^ 2;
t261 = t238 - t239;
t108 = t112 * qJD(3);
t260 = t108 * pkin(7) + t235;
t259 = -t105 * pkin(7) + t234;
t237 = t52 * t183;
t148 = sin(pkin(8));
t150 = cos(pkin(8));
t186 = g(1) * t150 + g(2) * t148;
t177 = t186 * t157;
t154 = sin(qJ(2));
t244 = g(3) * t154;
t168 = t177 + t244;
t211 = t154 * qJD(1);
t125 = qJD(2) * pkin(6) + t211;
t231 = qJD(2) * pkin(2);
t126 = -t210 - t231;
t145 = t153 ^ 2;
t146 = t156 ^ 2;
t218 = t145 + t146;
t192 = t218 * t157;
t258 = t125 * t192 + t126 * t154;
t140 = g(3) * t157;
t167 = t186 * t154 - t140;
t242 = t106 * pkin(7);
t190 = qJ(4) * qJD(2) + t125;
t95 = t190 * t156;
t73 = t147 * t95;
t230 = qJD(3) * pkin(3);
t94 = t190 * t153;
t77 = -t94 + t230;
t37 = t149 * t77 - t73;
t23 = qJD(3) * pkin(4) - t242 + t37;
t243 = t103 * pkin(7);
t229 = t149 * t95;
t38 = t147 * t77 + t229;
t25 = t38 - t243;
t209 = qJD(1) * qJD(2);
t111 = qJDD(2) * pkin(6) + t154 * qJDD(1) + t157 * t209;
t173 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + t111;
t180 = qJD(3) * t190;
t33 = qJDD(3) * pkin(3) - t173 * t153 - t156 * t180;
t34 = -t153 * t180 + t173 * t156;
t16 = -t147 * t34 + t149 * t33;
t6 = qJDD(3) * pkin(4) - t58 * pkin(7) + t16;
t17 = t147 * t33 + t149 * t34;
t7 = -t57 * pkin(7) + t17;
t1 = (qJD(5) * t23 + t7) * t155 + t152 * t6 - t25 * t214;
t144 = qJ(3) + pkin(9);
t139 = qJ(5) + t144;
t132 = sin(t139);
t133 = cos(t139);
t221 = t150 * t157;
t223 = t148 * t157;
t240 = t156 * pkin(3);
t136 = pkin(2) + t240;
t109 = -t136 * qJD(2) + qJD(4) - t210;
t61 = t103 * pkin(4) + t109;
t257 = -t61 * t52 - g(1) * (-t148 * t132 - t133 * t221) - g(2) * (t150 * t132 - t133 * t223) + t133 * t244 - t1;
t9 = t152 * t23 + t155 * t25;
t2 = -t9 * qJD(5) - t152 * t7 + t155 * t6;
t256 = t183 * t61 - g(1) * (-t132 * t221 + t148 * t133) - g(2) * (-t132 * t223 - t150 * t133) + t132 * t244 + t2;
t253 = t153 * t230 - t211;
t204 = t157 * qJDD(1);
t188 = t154 * t209 - t204;
t226 = qJDD(2) * pkin(2);
t110 = t188 - t226;
t158 = qJD(3) ^ 2;
t252 = -pkin(6) * t158 + t154 * (t186 + t209) - t110 - t140 + t226;
t251 = t106 ^ 2;
t250 = t57 * pkin(4);
t122 = t236 * t153;
t123 = t236 * t156;
t63 = -t149 * t122 - t147 * t123;
t47 = -t113 * pkin(7) + t63;
t64 = -t147 * t122 + t149 * t123;
t48 = -t112 * pkin(7) + t64;
t18 = -t152 * t48 + t155 * t47;
t249 = t18 * qJD(5) + t260 * t152 + t259 * t155;
t19 = t152 * t47 + t155 * t48;
t248 = -t19 * qJD(5) - t259 * t152 + t260 * t155;
t247 = pkin(3) * t147;
t241 = t153 * pkin(3);
t40 = -t149 * t94 - t73;
t39 = t147 * t94 - t229;
t26 = t39 + t243;
t27 = t40 - t242;
t134 = t149 * pkin(3) + pkin(4);
t97 = t155 * t134 - t152 * t247;
t233 = t97 * qJD(5) - t152 * t26 - t155 * t27;
t98 = t152 * t134 + t155 * t247;
t232 = -t98 * qJD(5) + t152 * t27 - t155 * t26;
t225 = t106 * t103;
t220 = qJDD(1) - g(3);
t219 = t145 - t146;
t159 = qJD(2) ^ 2;
t217 = t158 + t159;
t215 = qJD(2) * t154;
t212 = t109 * qJD(2);
t207 = qJDD(3) * t153;
t203 = t157 * qJDD(2);
t201 = t153 * t159 * t156;
t138 = cos(t144);
t119 = pkin(4) * t138 + t240;
t193 = t218 * t111;
t191 = t105 * pkin(4) + t253;
t189 = t218 * qJDD(2);
t187 = t153 * t197;
t185 = g(1) * t148 - g(2) * t150;
t92 = t113 * t154;
t93 = t112 * t154;
t41 = t152 * t93 - t155 * t92;
t42 = -t152 * t92 - t155 * t93;
t60 = -t152 * t112 + t155 * t113;
t176 = t185 * t156;
t164 = -pkin(6) * qJDD(3) + (t126 + t210 - t231) * qJD(3);
t62 = pkin(3) * t198 - t136 * qJDD(2) + qJDD(4) + t188;
t163 = -t126 * qJD(2) - t111 + t168;
t160 = t62 - t167;
t142 = -pkin(7) - t236;
t141 = qJDD(3) + qJDD(5);
t137 = sin(t144);
t118 = -pkin(4) * t137 - t241;
t117 = pkin(2) + t119;
t101 = t103 ^ 2;
t81 = t112 * pkin(4) - t136;
t67 = pkin(3) * t216 + t106 * pkin(4);
t59 = t155 * t112 + t152 * t113;
t44 = -qJD(2) * t174 - qJD(3) * t92;
t43 = -t157 * t106 + t154 * t108;
t24 = t62 + t250;
t21 = qJD(5) * t60 + t155 * t105 - t152 * t108;
t20 = t152 * t105 + t155 * t108 + t112 * t213 + t113 * t214;
t11 = -qJD(5) * t42 - t152 * t44 + t155 * t43;
t10 = qJD(5) * t41 + t152 * t43 + t155 * t44;
t8 = -t152 * t25 + t155 * t23;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t220, 0, 0, 0, 0, 0, 0, -t159 * t154 + t203, -qJDD(2) * t154 - t159 * t157, 0, -g(3) + (t154 ^ 2 + t157 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, (-0.2e1 * t198 + t205) * t157 + (-t217 * t156 - t207) * t154, (-qJDD(3) * t154 - 0.2e1 * t157 * t208) * t156 + (t217 * t154 - t203) * t153, t154 * t189 + t159 * t192, t258 * qJD(2) - t110 * t157 + t154 * t193 - g(3), 0, 0, 0, 0, 0, 0, t43 * qJD(3) - t92 * qJDD(3) + t103 * t215 - t157 * t57, -t44 * qJD(3) + t93 * qJDD(3) + t106 * t215 - t157 * t58, -t44 * t103 - t43 * t106 + t93 * t57 + t92 * t58, t154 * t212 - t62 * t157 - t16 * t92 - t17 * t93 + t37 * t43 + t38 * t44 - g(3), 0, 0, 0, 0, 0, 0, t11 * t143 + t41 * t141 + t157 * t165 - t215 * t52, -t10 * t143 - t42 * t141 - t157 * t175 - t183 * t215, t10 * t52 + t11 * t183 + t165 * t42 - t175 * t41, t1 * t42 + t9 * t10 + t8 * t11 - t24 * t157 + t2 * t41 + t61 * t215 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t204 + t167, -t220 * t154 + t177, 0, 0, t145 * qJDD(2) + 0.2e1 * t187, 0.2e1 * t153 * t205 - 0.2e1 * t219 * t208, t158 * t156 + t207, t146 * qJDD(2) - 0.2e1 * t187, qJDD(3) * t156 - t158 * t153, 0, t164 * t153 + t252 * t156, -t252 * t153 + t164 * t156, -t244 + t193 + pkin(6) * t189 + (-t218 * t209 - t186) * t157, (-t110 + t167) * pkin(2) + (t193 - t168) * pkin(6) - t258 * qJD(1), -t106 * t108 + t58 * t113, t108 * t103 - t106 * t105 - t58 * t112 - t113 * t57, -t108 * qJD(3) + t113 * qJDD(3), t103 * t105 + t57 * t112, -t105 * qJD(3) - t112 * qJDD(3), 0, -t103 * t211 + t63 * qJDD(3) + t109 * t105 + t62 * t112 - t136 * t57 + t167 * t138 + (t103 * t241 + t235) * qJD(3), -t106 * t211 - t64 * qJDD(3) - t109 * t108 + t62 * t113 - t136 * t58 - t167 * t137 + (t106 * t241 - t234) * qJD(3), -t234 * t103 - t38 * t105 - t235 * t106 + t37 * t108 - t17 * t112 - t16 * t113 - t64 * t57 - t63 * t58 - t168, t17 * t64 + t16 * t63 - t62 * t136 - g(3) * (t157 * t136 + t154 * t236) + t234 * t38 + t235 * t37 + t253 * t109 + t186 * (t136 * t154 - t157 * t236), t175 * t60 + t183 * t20, t165 * t60 - t175 * t59 + t183 * t21 - t20 * t52, t60 * t141 - t20 * t143, -t165 * t59 - t21 * t52, -t59 * t141 - t21 * t143, 0, t167 * t133 + t18 * t141 + t248 * t143 - t165 * t81 - t191 * t52 + t61 * t21 + t24 * t59, -t132 * t167 - t19 * t141 - t249 * t143 + t175 * t81 - t183 * t191 - t61 * t20 + t24 * t60, -t1 * t59 + t165 * t19 - t175 * t18 + t183 * t248 - t2 * t60 + t8 * t20 - t9 * t21 + t249 * t52 - t168, t1 * t19 + t2 * t18 + t24 * t81 - g(3) * (t157 * t117 - t154 * t142) + t249 * t9 + t248 * t8 + t191 * t61 + t186 * (t117 * t154 + t142 * t157); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201, t219 * t159, t206, t201, t205, qJDD(3), t163 * t153 - t176, t185 * t153 + t156 * t163, 0, 0, t225, -t101 + t251, (t103 + t200) * qJD(3) + t170, -t225, -t184, qJDD(3), -t39 * qJD(3) - t109 * t106 - g(1) * (-t137 * t221 + t148 * t138) - g(2) * (-t137 * t223 - t150 * t138) + t137 * t244 + (t149 * qJDD(3) - t103 * t216) * pkin(3) + t16, t40 * qJD(3) + t109 * t103 - g(1) * (-t148 * t137 - t138 * t221) - g(2) * (t150 * t137 - t138 * t223) + t138 * t244 + (-t147 * qJDD(3) - t106 * t216) * pkin(3) - t17, (t38 + t39) * t106 + (-t37 + t40) * t103 + (-t147 * t57 - t149 * t58) * pkin(3), -t37 * t39 - t38 * t40 + (t17 * t147 + t16 * t149 - t176 + (t168 - t212) * t153) * pkin(3), t237, t261, t262, -t237, t263, t141, t97 * t141 + t232 * t143 + t52 * t67 + t256, -t98 * t141 - t233 * t143 + t183 * t67 + t257, t98 * t165 - t175 * t97 + (t233 + t8) * t52 + (t232 - t9) * t183, t1 * t98 + t2 * t97 - t61 * t67 - g(1) * (t118 * t221 + t148 * t119) - g(2) * (t118 * t223 - t150 * t119) - t118 * t244 + t233 * t9 + t232 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t106 * qJD(3) + t184, (-t103 + t200) * qJD(3) + t170, -t101 - t251, t38 * t103 + t37 * t106 + t160, 0, 0, 0, 0, 0, 0, -t165 - t228, t175 + t227, -t238 - t239, -t183 * t8 - t9 * t52 + t160 + t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, t261, t262, -t237, t263, t141, t9 * t143 + t256, t8 * t143 + t257, 0, 0;];
tau_reg = t3;
