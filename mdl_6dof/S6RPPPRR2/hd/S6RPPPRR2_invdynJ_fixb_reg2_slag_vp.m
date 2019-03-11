% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPPRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:08
% EndTime: 2019-03-09 01:32:12
% DurationCPUTime: 2.56s
% Computational Cost: add. (4686->351), mult. (8630->429), div. (0->0), fcn. (5917->14), ass. (0->196)
t132 = sin(pkin(10));
t134 = cos(pkin(10));
t138 = sin(qJ(5));
t141 = cos(qJ(5));
t86 = t132 * t141 + t134 * t138;
t78 = t86 * qJD(1);
t253 = qJD(6) + t78;
t137 = sin(qJ(6));
t140 = cos(qJ(6));
t208 = qJD(1) * t134;
t135 = cos(pkin(9));
t107 = -pkin(1) * t135 - pkin(2);
t99 = -qJ(4) + t107;
t85 = qJD(1) * t99 + qJD(3);
t60 = -qJD(2) * t132 + t134 * t85;
t52 = -pkin(7) * t208 + t60;
t209 = qJD(1) * t132;
t61 = qJD(2) * t134 + t132 * t85;
t53 = -pkin(7) * t209 + t61;
t25 = t138 * t52 + t141 * t53;
t23 = qJD(5) * pkin(8) + t25;
t133 = sin(pkin(9));
t102 = pkin(1) * t133 + qJ(3);
t90 = t102 * qJD(1);
t88 = qJD(4) + t90;
t74 = pkin(4) * t209 + t88;
t192 = t141 * t208;
t193 = t138 * t209;
t80 = t192 - t193;
t28 = pkin(5) * t78 - pkin(8) * t80 + t74;
t9 = -t137 * t23 + t140 * t28;
t255 = t9 * t253;
t10 = t137 * t28 + t140 * t23;
t254 = t10 * t253;
t126 = pkin(10) + qJ(5);
t115 = sin(t126);
t117 = cos(t126);
t127 = qJ(1) + pkin(9);
t116 = sin(t127);
t118 = cos(t127);
t247 = g(1) * t116 - g(2) * t118;
t151 = g(3) * t115 - t117 * t247;
t200 = t134 * qJDD(1);
t202 = qJD(1) * qJD(4);
t70 = qJDD(1) * t99 + qJDD(3) - t202;
t55 = -qJDD(2) * t132 + t134 * t70;
t49 = -pkin(7) * t200 + t55;
t201 = t132 * qJDD(1);
t56 = qJDD(2) * t134 + t132 * t70;
t50 = -pkin(7) * t201 + t56;
t190 = t138 * t50 - t141 * t49;
t8 = -t25 * qJD(5) - t190;
t6 = -qJDD(5) * pkin(5) - t8;
t149 = -t6 + t151;
t252 = -pkin(8) * qJD(6) * t253 + t149;
t186 = t140 * t253;
t153 = -qJD(5) * t193 + qJDD(1) * t86;
t45 = qJD(5) * t192 + t153;
t42 = qJDD(6) + t45;
t225 = t137 * t42;
t251 = -t186 * t253 - t225;
t203 = t140 * qJD(5);
t205 = qJD(6) * t137;
t167 = -t138 * t201 + t141 * t200;
t82 = t86 * qJD(5);
t44 = qJD(1) * t82 - t167;
t20 = -qJD(6) * t203 - qJDD(5) * t137 + t140 * t44 + t205 * t80;
t59 = qJD(5) * t137 + t140 * t80;
t206 = qJD(5) * t141;
t207 = qJD(5) * t138;
t83 = -t132 * t207 + t134 * t206;
t228 = -t20 * t86 + t59 * t83;
t187 = t137 * t253;
t249 = t59 * t187;
t222 = pkin(1) * qJDD(1);
t177 = g(1) * t118 + g(2) * t116;
t248 = t177 * t117;
t220 = qJD(6) * t59;
t21 = -qJDD(5) * t140 - t137 * t44 + t220;
t87 = -t132 * t138 + t134 * t141;
t246 = t107 * qJDD(1);
t108 = pkin(4) * t201;
t129 = qJD(3) * qJD(1);
t109 = t133 * t222;
t211 = qJDD(1) * qJ(3) + t109;
t194 = -t129 - t211;
t84 = qJDD(4) - t194;
t68 = t108 + t84;
t14 = pkin(5) * t45 + pkin(8) * t44 + t68;
t13 = t140 * t14;
t199 = -t138 * t49 - t141 * t50 - t206 * t52;
t7 = -t207 * t53 - t199;
t5 = qJDD(5) * pkin(8) + t7;
t2 = -qJD(6) * t10 - t137 * t5 + t13;
t245 = t2 + t254;
t240 = g(3) * t117;
t150 = -t115 * t247 - t240;
t171 = -t253 * t82 + t42 * t87;
t197 = t87 * t205;
t244 = -t140 * t171 + t197 * t253;
t243 = t80 ^ 2;
t242 = -pkin(7) + t99;
t239 = t10 * t83;
t122 = t132 * pkin(4);
t57 = t137 * t80 - t203;
t237 = t57 * t78;
t236 = t59 * t57;
t235 = t59 * t80;
t234 = t80 * t57;
t233 = t80 * t78;
t232 = t82 * t57;
t231 = t82 * t59;
t204 = qJD(6) * t140;
t230 = -t137 * t21 - t204 * t57;
t223 = t21 * t140;
t229 = t140 * t232 - t223 * t87;
t227 = -t45 * t87 + t78 * t82;
t226 = qJD(6) * t9;
t224 = t138 * t53;
t35 = t140 * t42;
t221 = qJD(6) * t57;
t219 = t116 * t137;
t218 = t116 * t140;
t217 = t118 * t137;
t216 = t118 * t140;
t124 = t132 ^ 2;
t125 = t134 ^ 2;
t210 = t124 + t125;
t196 = t87 * t204;
t142 = cos(qJ(1));
t195 = pkin(1) * t142 + pkin(2) * t118 + qJ(3) * t116;
t139 = sin(qJ(1));
t191 = -t139 * pkin(1) + qJ(3) * t118;
t189 = -t20 + t221;
t185 = qJD(6) * t86 + qJD(1);
t184 = t210 * qJDD(1);
t182 = t59 * t196;
t181 = -t2 * t86 - t9 * t83;
t89 = t102 + t122;
t24 = t141 * t52 - t224;
t22 = -qJD(5) * pkin(5) - t24;
t180 = -t22 * t82 + t6 * t87;
t1 = t137 * t14 + t140 * t5 + t226;
t179 = t1 - t255;
t178 = pkin(5) * t115 - pkin(8) * t117;
t175 = g(1) * t139 - g(2) * t142;
t174 = -t87 * t20 - t231;
t173 = t87 * t21 - t232;
t172 = -t86 * t21 - t83 * t57;
t170 = -t44 * t87 - t80 * t82;
t169 = -t44 * t86 + t80 * t83;
t168 = t45 * t86 + t78 * t83;
t166 = -t10 * t140 + t137 * t9;
t165 = t10 * t137 + t140 * t9;
t164 = t56 * t132 + t55 * t134;
t163 = t132 * t61 + t134 * t60;
t32 = pkin(5) * t86 - pkin(8) * t87 + t89;
t76 = t242 * t132;
t77 = t242 * t134;
t34 = t138 * t77 + t141 * t76;
t15 = -t137 * t34 + t140 * t32;
t16 = t137 * t32 + t140 * t34;
t33 = t138 * t76 - t141 * t77;
t46 = -qJD(5) * t82 + qJDD(5) * t87;
t47 = -qJD(5) * t83 - qJDD(5) * t86;
t162 = t35 + (-t137 * t78 - t205) * t253;
t161 = -t116 * pkin(2) + t191;
t160 = -qJD(6) * t28 + t240 - t5;
t159 = -pkin(8) * t42 + t22 * t253;
t136 = -pkin(7) - qJ(4);
t158 = t116 * t122 - t118 * t136 + t195;
t157 = qJDD(3) + t246;
t156 = -t164 + t247;
t155 = t116 * t136 + t118 * t122 + t161;
t154 = qJDD(1) * t102 - t177;
t152 = -t177 + t84;
t148 = -t24 * t82 + t25 * t83 + t7 * t86 + t8 * t87 - t247;
t147 = -t137 * t171 - t196 * t253;
t146 = t154 + t84 + t129;
t145 = -qJD(6) * t165 + t1 * t140 - t2 * t137;
t143 = qJD(1) ^ 2;
t131 = qJDD(2) - g(3);
t75 = t78 ^ 2;
t67 = t115 * t216 - t219;
t66 = t115 * t217 + t218;
t65 = t115 * t218 + t217;
t64 = -t115 * t219 + t216;
t43 = pkin(5) * t80 + pkin(8) * t78;
t40 = pkin(5) * t83 + pkin(8) * t82 + qJD(3);
t27 = qJD(4) * t87 + qJD(5) * t34;
t26 = -qJD(4) * t86 - qJD(5) * t33;
t12 = t137 * t43 + t140 * t24;
t11 = -t137 * t24 + t140 * t43;
t4 = -qJD(6) * t16 - t137 * t26 + t140 * t40;
t3 = qJD(6) * t15 + t137 * t40 + t140 * t26;
t17 = [0, 0, 0, 0, 0, qJDD(1), t175, g(1) * t142 + g(2) * t139, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t135 * t222 + t247, -0.2e1 * t109 + t177, 0 (t175 + (t133 ^ 2 + t135 ^ 2) * t222) * pkin(1), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - t247 + 0.2e1 * t246, 0.2e1 * t129 + t154 + t211, -g(1) * t161 - g(2) * t195 + t90 * qJD(3) - t102 * t194 + t107 * t157, t125 * qJDD(1), -0.2e1 * t132 * t200, 0, t124 * qJDD(1), 0, 0, t146 * t132, t146 * t134, -t184 * t99 + t202 * t210 + t156, t84 * t102 + t88 * qJD(3) - g(1) * ((-pkin(2) - qJ(4)) * t116 + t191) - g(2) * (qJ(4) * t118 + t195) + t164 * t99 - t163 * qJD(4), t170, -t169 + t227, t46, t168, t47, 0, qJD(3) * t78 - t27 * qJD(5) - t33 * qJDD(5) - t115 * t177 + t89 * t45 + t68 * t86 + t74 * t83, qJD(3) * t80 - t26 * qJD(5) - t34 * qJDD(5) - t89 * t44 + t68 * t87 - t74 * t82 - t248, -t26 * t78 + t27 * t80 - t33 * t44 - t34 * t45 - t148, -g(1) * t155 - g(2) * t158 + qJD(3) * t74 - t24 * t27 + t25 * t26 - t33 * t8 + t34 * t7 + t68 * t89, t140 * t174 - t197 * t59, -t182 + (t231 + (t20 + t221) * t87) * t137 + t229, t228 - t244, t137 * t173 + t196 * t57, t147 + t172, t253 * t83 + t42 * t86, -g(1) * t67 - g(2) * t65 + t137 * t180 + t15 * t42 + t196 * t22 + t33 * t21 + t253 * t4 + t27 * t57 - t181, g(1) * t66 - g(2) * t64 - t1 * t86 + t140 * t180 - t16 * t42 - t197 * t22 - t33 * t20 - t253 * t3 + t27 * t59 - t239, t15 * t20 - t16 * t21 - t3 * t57 - t4 * t59 + t165 * t82 + t248 + (qJD(6) * t166 - t1 * t137 - t2 * t140) * t87, t1 * t16 + t10 * t3 + t2 * t15 + t9 * t4 + t6 * t33 + t22 * t27 - g(1) * (t118 * t178 + t155) - g(2) * (t116 * t178 + t158); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132 * t55 + t134 * t56 - g(3), 0, 0, 0, 0, 0, 0, t47, -t46, t169 + t227, -t24 * t83 - t25 * t82 + t7 * t87 - t8 * t86 - g(3), 0, 0, 0, 0, 0, 0, t147 - t172, t228 + t244, t182 + (t189 * t87 - t231) * t137 + t229, t145 * t87 + t166 * t82 + t22 * t83 + t6 * t86 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t143, -t90 * qJD(1) + t157 - t247, 0, 0, 0, 0, 0, 0, -t143 * t132, -t143 * t134, -t184, -qJD(1) * t88 - t156, 0, 0, 0, 0, 0, 0, -qJD(1) * t78 + t46, -qJD(1) * t80 + t47, -t168 - t170, -qJD(1) * t74 + t148, 0, 0, 0, 0, 0, 0, -t86 * t225 + (-t137 * t83 - t140 * t185) * t253 - t173, -t86 * t35 + (t137 * t185 - t140 * t83) * t253 - t174 (t185 * t59 + t172) * t140 + (t185 * t57 + t228) * t137 (-qJD(1) * t9 + t239 + (t1 - t226) * t86) * t140 + (-t10 * t185 + t181) * t137 - t180 - t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t200, -t210 * t143, qJD(1) * t163 + t152, 0, 0, 0, 0, 0, 0 (t80 + t192) * qJD(5) + t153, -0.2e1 * t78 * qJD(5) + t167, -t75 - t243, t24 * t80 + t25 * t78 + t108 + t152, 0, 0, 0, 0, 0, 0, t162 - t234, -t235 + t251 (t20 - t237) * t140 + t249 + t230, t137 * t179 + t140 * t245 - t22 * t80 - t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, -t75 + t243, t167, -t233 (t80 - t192) * qJD(5) - t153, qJDD(5), -t74 * t80 + t151 - t190, t74 * t78 + (t24 + t224) * qJD(5) + t199 - t150, 0, 0, -t137 * t20 + t186 * t59 (-t20 - t237) * t140 - t249 + t230, -t235 - t251, t187 * t57 - t223, t162 + t234, -t253 * t80, -pkin(5) * t21 - t11 * t253 + t137 * t159 + t140 * t252 - t25 * t57 - t9 * t80, pkin(5) * t20 + t10 * t80 + t12 * t253 - t137 * t252 + t140 * t159 - t25 * t59, t11 * t59 + t12 * t57 + ((-t21 + t220) * pkin(8) + t179) * t140 + (pkin(8) * t189 - t245) * t137 + t150, -t10 * t12 - t9 * t11 - t22 * t25 + t149 * pkin(5) + (t145 + t150) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, -t57 ^ 2 + t59 ^ 2, t253 * t57 - t20, -t236, t253 * t59 - t21, t42, -g(1) * t64 - g(2) * t66 + t137 * t160 - t204 * t23 - t22 * t59 + t13 + t254, g(1) * t65 - g(2) * t67 + t22 * t57 + t255 + (qJD(6) * t23 - t14) * t137 + t160 * t140, 0, 0;];
tau_reg  = t17;
