% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:46
% EndTime: 2019-12-31 20:53:50
% DurationCPUTime: 2.18s
% Computational Cost: add. (1837->345), mult. (2682->355), div. (0->0), fcn. (1322->8), ass. (0->192)
t127 = qJDD(1) + qJDD(2);
t135 = sin(qJ(3));
t118 = t135 * qJ(4);
t138 = cos(qJ(3));
t208 = t138 * pkin(3) + t118;
t255 = -pkin(2) - t208;
t257 = t127 * t255;
t128 = qJD(1) + qJD(2);
t256 = t128 * t255;
t131 = t135 ^ 2;
t132 = t138 ^ 2;
t207 = t131 + t132;
t136 = sin(qJ(2));
t198 = qJDD(1) * t136;
t139 = cos(qJ(2));
t204 = qJD(2) * t139;
t40 = t127 * pkin(7) + (qJD(1) * t204 + t198) * pkin(1);
t254 = t207 * t40;
t213 = t128 * t138;
t243 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t215 = t127 * t132;
t216 = t127 * t131;
t253 = t215 + t216;
t182 = t207 * t139;
t228 = pkin(1) * qJD(1);
t194 = t136 * t228;
t65 = pkin(7) * t128 + t194;
t206 = qJD(1) * t139;
t192 = pkin(1) * t206;
t234 = t128 * pkin(2);
t66 = -t192 - t234;
t252 = t136 * t66 + t182 * t65;
t232 = pkin(3) + qJ(5);
t251 = t232 * qJD(3);
t250 = t232 * qJDD(3);
t200 = qJD(3) * t138;
t188 = t128 * t200;
t96 = t135 * t127;
t248 = (-t188 - t96) * pkin(4);
t33 = t135 * t40;
t46 = t65 * t200;
t191 = qJDD(4) + t33 + t46;
t223 = qJDD(3) * pkin(3);
t10 = t191 - t223;
t201 = qJD(3) * t135;
t35 = t138 * t40;
t9 = t65 * t201 - t243 - t35;
t247 = t10 * t135 - t9 * t138;
t203 = qJD(3) * qJ(4);
t177 = -qJD(5) - t203;
t51 = t138 * t65;
t30 = pkin(4) * t213 + t51;
t20 = -t177 + t30;
t50 = t135 * t65;
t245 = -qJD(4) - t50;
t133 = qJ(1) + qJ(2);
t116 = sin(t133);
t117 = cos(t133);
t209 = g(1) * t117 + g(2) * t116;
t244 = 0.2e1 * t243;
t219 = t117 * t135;
t221 = t116 * t135;
t242 = g(1) * t219 + g(2) * t221 - g(3) * t138;
t37 = -qJD(3) * pkin(3) - t245;
t43 = -t51 - t203;
t241 = t37 * t200 + t43 * t201 + t247;
t202 = qJD(3) * t128;
t97 = t138 * t127;
t240 = 0.2e1 * t135 * t97 + 0.2e1 * (-t131 + t132) * t202;
t239 = pkin(1) * t139;
t141 = qJD(3) ^ 2;
t238 = pkin(7) * t141;
t106 = g(1) * t116;
t137 = sin(qJ(1));
t237 = g(1) * t137;
t236 = g(3) * t135;
t235 = t127 * pkin(2);
t231 = -g(1) * t221 + g(2) * t219;
t218 = t117 * t138;
t220 = t116 * t138;
t230 = g(1) * t220 - g(2) * t218;
t229 = qJ(4) * t97 + qJD(4) * t213;
t226 = t37 * t138;
t205 = qJD(2) * t136;
t114 = pkin(1) * t205;
t225 = qJD(1) * t114 - qJDD(1) * t239;
t224 = pkin(7) * qJDD(3);
t109 = pkin(1) * t136 + pkin(7);
t222 = t109 * t141;
t126 = t128 ^ 2;
t217 = t126 * t132;
t214 = t128 * t135;
t29 = -pkin(4) * t214 - t50;
t212 = qJD(4) - t29;
t102 = t117 * pkin(7);
t211 = t117 * pkin(4) + t102;
t210 = t117 * pkin(2) + t116 * pkin(7);
t199 = t135 * qJD(4);
t197 = qJDD(3) * t109;
t179 = t128 * t194;
t62 = t135 * t179;
t196 = t62 - t231;
t112 = pkin(3) * t201;
t195 = t128 * t112 + t225;
t193 = pkin(1) * t204;
t190 = pkin(4) * t97 + qJDD(5) + t35;
t110 = -pkin(2) - t239;
t189 = t128 * t205;
t187 = -pkin(4) * t128 - t65;
t22 = -t192 + t256;
t186 = -t22 - t256;
t39 = t225 - t235;
t184 = t39 * t135 + t66 * t200 + t231;
t183 = -t138 * t179 - t192 * t201 - t230;
t181 = -t33 + t242;
t180 = pkin(3) * t218 + t117 * t118 + t210;
t176 = t135 * t188;
t175 = g(1) * (-t116 * pkin(2) + t102);
t174 = -t235 + t238;
t140 = cos(qJ(1));
t172 = -g(2) * t140 + t237;
t171 = qJ(5) * t138 + t208;
t170 = t187 * qJD(3);
t169 = -qJDD(4) + t181;
t167 = t43 * t135 + t226;
t166 = t135 * t37 - t138 * t43;
t165 = -qJ(4) * t138 + qJ(5) * t135;
t164 = -g(2) * t117 + t106 - t225;
t49 = -pkin(2) - t171;
t163 = t182 * t228;
t162 = t207 * t128 * t193 + t253 * t109 - t209;
t161 = t116 * pkin(4) + qJ(5) * t218 + t180;
t155 = -t232 * t138 - pkin(2) - t118;
t1 = (t165 * qJD(3) - qJD(5) * t138 - t199) * t128 + t155 * t127 + t195;
t21 = qJ(5) * t201 + t177 * t138 + t112 - t199;
t16 = t114 + t21;
t38 = t49 - t239;
t159 = -t127 * t38 - t128 * t16 - t1;
t158 = -t127 * t49 - t128 * t21 - t1;
t157 = -t46 + t169;
t156 = -qJ(4) * t200 - t199;
t44 = t112 + t156;
t154 = t128 * t44 + t238 + t257;
t31 = t114 + t44;
t52 = t110 - t208;
t153 = t127 * t52 + t128 * t31 + t222;
t17 = t212 - t251;
t151 = t191 - t248 - t250;
t4 = -qJD(3) * qJD(5) + t151;
t6 = t135 * t170 + t190 + t243;
t152 = t4 * t135 + t6 * t138 + t17 * t200 - t20 * t201 - t209;
t150 = t155 * t106;
t149 = pkin(1) * t189 + t110 * t127 + t222;
t14 = t155 * t128 - t192;
t148 = -t138 * t209 + t14 * t213 + t190;
t147 = -g(1) * t102 - t106 * t255;
t146 = t253 * pkin(7) - t128 * t163 - t209;
t145 = -t197 + (t110 * t128 - t193) * qJD(3);
t144 = t197 + (-t128 * t52 + t193 - t22) * qJD(3);
t143 = t167 * qJD(3) + t247;
t124 = t140 * pkin(1);
t122 = t138 * pkin(4);
t121 = t135 * pkin(4);
t113 = pkin(4) * t200;
t99 = t131 * t126;
t95 = pkin(3) * t214;
t80 = qJ(4) * t218;
t77 = qJ(4) * t220;
t74 = pkin(7) * t138 + t122;
t73 = pkin(7) * t135 + t121;
t72 = t135 * t126 * t138;
t70 = -t99 - t141;
t69 = qJDD(3) * t138 - t141 * t135;
t68 = qJDD(3) * t135 + t138 * t141;
t60 = qJDD(3) + t72;
t59 = pkin(7) * t200 + t113;
t58 = (-pkin(4) - pkin(7)) * t201;
t55 = t109 * t138 + t122;
t54 = t109 * t135 + t121;
t53 = -t99 + t217;
t47 = t66 * t201;
t45 = -qJ(4) * t213 + t95;
t42 = -0.2e1 * t176 + t215;
t41 = 0.2e1 * t176 + t216;
t26 = t165 * t128 + t95;
t25 = t109 * t200 + t135 * t193 + t113;
t24 = t138 * t193 + (-pkin(4) - t109) * t201;
t15 = t22 * t214;
t12 = t14 * t201;
t7 = t156 * t128 + t195 + t257;
t5 = t7 * t138;
t2 = [0, 0, 0, 0, 0, qJDD(1), t172, g(1) * t140 + g(2) * t137, 0, 0, 0, 0, 0, 0, 0, t127, (t127 * t139 - t189) * pkin(1) + t164, ((-qJDD(1) - t127) * t136 + (-qJD(1) - t128) * t204) * pkin(1) + t209, 0, (t172 + (t136 ^ 2 + t139 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t41, t240, t68, t42, t69, 0, t47 + t145 * t135 + (-t149 - t39) * t138 + t230, t135 * t149 + t138 * t145 + t184, t162 + t254, t39 * t110 - t175 - g(2) * (t124 + t210) + t109 * t254 + (t252 * qJD(2) + t237) * pkin(1), 0, -t68, -t69, t41, t240, t42, t162 + t241, t135 * t144 + t138 * t153 - t230 + t5, t144 * t138 + (-t153 - t7) * t135 - t231, t7 * t52 + t22 * t31 - g(2) * (t124 + t180) + (t166 * t204 + t237) * pkin(1) + t143 * t109 + t147, 0, -t69, t68, t42, -t240, t41, (t135 * t54 + t138 * t55) * t127 + (t135 * t25 + t138 * t24 + (-t135 * t55 + t138 * t54) * qJD(3)) * t128 + t152, t55 * qJDD(3) + t159 * t135 + (t24 + (-t128 * t38 - t14) * t138) * qJD(3) - t231, -t54 * qJDD(3) + t12 + (t38 * t214 - t25) * qJD(3) + t159 * t138 + t230, t1 * t38 + t14 * t16 + t4 * t54 + t17 * t25 + t6 * t55 + t20 * t24 - g(1) * (-t137 * pkin(1) + t211) - g(2) * (t124 + t161) - t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t164 + t179, (-t198 + (-qJD(2) + t128) * t206) * pkin(1) + t209, 0, 0, t41, t240, t68, t42, t69, 0, t47 + (-pkin(2) * t202 - t224) * t135 + (-t174 - t39) * t138 - t183, -t62 + t174 * t135 + (-t224 + (t192 - t234) * qJD(3)) * t138 + t184, t146 + t254, -t39 * pkin(2) + pkin(7) * t254 - g(2) * t210 - t252 * t228 - t175, 0, -t68, -t69, t41, t240, t42, t146 + t241, t5 + t154 * t138 + (qJD(3) * t186 + t224) * t135 + t183, (t224 + (t186 - t192) * qJD(3)) * t138 + (-t154 - t7) * t135 + t196, t7 * t255 + t22 * t44 - g(2) * t180 + (-t136 * t22 - t139 * t166) * t228 + t143 * pkin(7) + t147, 0, -t69, t68, t42, -t240, t41, (t135 * t73 + t138 * t74) * t127 + (t135 * t59 + t138 * t58 + (-t135 * t74 + t138 * t73) * qJD(3) - t163) * t128 + t152, qJDD(3) * t74 + t158 * t135 + (t58 + (-t128 * t49 - t14 - t192) * t138) * qJD(3) + t196, -qJDD(3) * t73 + t12 + (t49 * t214 - t59) * qJD(3) + t158 * t138 - t183, t1 * t49 + t14 * t21 + t4 * t73 + t17 * t59 + t6 * t74 + t20 * t58 - g(1) * t211 - g(2) * t161 - t150 + (-t136 * t14 + (-t135 * t17 - t138 * t20) * t139) * t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t53, t96, t72, t97, qJDD(3), -t66 * t214 + t181, t236 - t35 + (-t128 * t66 + t209) * t138, 0, 0, qJDD(3), -t96, -t97, -t72, -t53, t72, -pkin(3) * t96 + (-qJD(3) * t208 - t167) * t128 + t229, -t45 * t213 + t15 - t169 - 0.2e1 * t223, t35 + (t128 * t45 - g(3)) * t135 + (t128 * t22 - t209) * t138 + t244, -t9 * qJ(4) - t10 * pkin(3) - t22 * t45 - t65 * t226 - g(1) * (-pkin(3) * t219 + t80) - g(2) * (-pkin(3) * t221 + t77) - g(3) * t208 + t245 * t43, qJDD(3), -t97, t96, t72, t53, -t72, -t232 * t96 + (-t17 - t29 - t251) * t213 + t229, -qJD(3) * t29 + (t128 * t26 - g(3) + t170) * t135 + t148 + t244, (-t135 * t14 + t138 * t26) * t128 + (0.2e1 * qJD(5) + t30) * qJD(3) + 0.2e1 * t250 + t157 + t248, t6 * qJ(4) - t14 * t26 - g(1) * t80 - g(2) * t77 - g(3) * t171 + t212 * t20 + (-qJD(5) - t30) * t17 + (t209 * t135 - t4) * t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t60, t70, qJD(3) * t43 + t15 - t157 - t223, 0, 0, 0, 0, 0, 0, t96, t70, -t60, t14 * t214 + (-qJD(5) - t20) * qJD(3) + t151 - t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, qJDD(3) - t72, -t141 - t217, -t236 + (t135 * t187 + t17) * qJD(3) + t148 + t243;];
tau_reg = t2;
