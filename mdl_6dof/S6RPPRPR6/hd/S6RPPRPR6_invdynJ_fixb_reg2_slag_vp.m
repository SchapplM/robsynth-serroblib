% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:40
% EndTime: 2019-03-09 01:51:44
% DurationCPUTime: 3.05s
% Computational Cost: add. (2512->400), mult. (4203->444), div. (0->0), fcn. (2119->6), ass. (0->212)
t122 = sin(qJ(1));
t125 = cos(qJ(1));
t208 = g(1) * t125 + g(2) * t122;
t119 = pkin(1) + qJ(3);
t256 = qJD(1) * t119;
t69 = -qJD(2) + t256;
t211 = t69 * qJD(1);
t263 = t208 + t211;
t120 = sin(qJ(6));
t123 = cos(qJ(6));
t199 = qJD(4) * t123;
t121 = sin(qJ(4));
t205 = qJD(1) * t121;
t50 = t120 * t205 + t199;
t222 = qJD(6) * t50;
t124 = cos(qJ(4));
t194 = qJD(1) * qJD(4);
t177 = t124 * t194;
t189 = t121 * qJDD(1);
t260 = t177 + t189;
t15 = t120 * qJDD(4) - t260 * t123 + t222;
t203 = qJD(1) * t124;
t79 = qJD(6) + t203;
t243 = t50 * t79;
t262 = t15 - t243;
t92 = qJ(2) * qJD(1) + qJD(3);
t68 = -pkin(7) * qJD(1) + t92;
t31 = (pkin(5) * qJD(1) - t68) * t124;
t210 = qJD(5) + t31;
t192 = qJDD(4) * qJ(5);
t112 = qJD(1) * qJD(2);
t113 = qJ(2) * qJDD(1);
t175 = qJDD(3) + t112 + t113;
t51 = -pkin(7) * qJDD(1) + t175;
t37 = t121 * t51;
t167 = -t37 - t192;
t231 = t124 * t68;
t16 = (-qJD(5) - t231) * qJD(4) + t167;
t166 = -t124 * t51 + qJDD(5);
t200 = qJD(4) * t121;
t47 = t68 * t200;
t151 = t166 + t47;
t221 = qJDD(4) * pkin(4);
t17 = t151 - t221;
t169 = -qJD(5) + t231;
t239 = qJD(4) * pkin(4);
t33 = -t169 - t239;
t228 = t33 * t121;
t202 = qJD(4) * qJ(5);
t56 = t121 * t68;
t36 = -t56 - t202;
t130 = (-t36 * t124 + t228) * qJD(4) - t16 * t121 - t17 * t124;
t261 = t130 - t208;
t247 = pkin(8) * t121;
t98 = t124 * qJ(5);
t148 = -t98 + t247;
t206 = pkin(4) * t205 - qJD(2);
t20 = (t148 + t119) * qJD(1) + t206;
t126 = -pkin(4) - pkin(8);
t21 = t126 * qJD(4) + t210;
t6 = t120 * t21 + t123 * t20;
t110 = qJDD(1) * qJ(3);
t111 = qJD(3) * qJD(1);
t117 = qJDD(1) * pkin(1);
t188 = t117 - qJDD(2);
t160 = -t110 - t111 - t188;
t178 = t121 * t194;
t134 = t260 * pkin(4) + qJ(5) * t178 - t160;
t170 = qJD(4) * pkin(8) - qJD(5);
t195 = qJ(5) * qJDD(1);
t7 = pkin(8) * t189 + (qJD(1) * t170 - t195) * t124 + t134;
t94 = t124 * qJDD(1);
t257 = -t178 + t94;
t9 = t257 * pkin(5) + t126 * qJDD(4) + t151;
t2 = -qJD(6) * t6 - t120 * t7 + t123 * t9;
t259 = t6 * t79 + t2;
t204 = qJD(1) * t123;
t179 = t121 * t204;
t197 = qJD(6) * t120;
t14 = qJD(4) * t197 - qJD(6) * t179 - t123 * qJDD(4) - t260 * t120;
t201 = qJD(4) * t120;
t48 = -t179 + t201;
t142 = t48 * t79;
t135 = t14 - t142;
t46 = -qJDD(6) - t257;
t34 = t123 * t46;
t138 = -t79 * t197 - t34;
t115 = t121 ^ 2;
t116 = t124 ^ 2;
t207 = t115 + t116;
t30 = -pkin(5) * t205 + t56;
t25 = t30 + t202;
t255 = -t126 * t46 + t25 * t79;
t118 = -pkin(7) + qJ(2);
t191 = qJDD(4) * t118;
t254 = qJD(4) * (qJD(2) + t69 + t256) + t191;
t26 = (-t98 + t119) * qJD(1) + t206;
t103 = t121 * pkin(4);
t226 = t103 - t98;
t57 = t226 + t119;
t253 = (qJD(1) * t57 + qJD(2) + t26) * qJD(4) + t191;
t5 = -t120 * t20 + t123 * t21;
t1 = qJD(6) * t5 + t120 * t9 + t123 * t7;
t155 = t120 * t5 - t123 * t6;
t252 = -qJD(6) * t155 + t1 * t120 + t2 * t123;
t156 = t120 * t6 + t123 * t5;
t251 = t156 * qJD(6) - t1 * t123 + t2 * t120;
t101 = 0.2e1 * t112;
t250 = pkin(5) + pkin(7);
t249 = t5 * t79;
t108 = g(1) * t122;
t107 = g(2) * t125;
t105 = g(3) * t121;
t246 = g(3) * t124;
t245 = t125 * pkin(7);
t244 = t50 * t48;
t242 = pkin(5) - t118;
t216 = t122 * t124;
t219 = t121 * t122;
t241 = pkin(4) * t216 + qJ(5) * t219;
t214 = t124 * t125;
t218 = t121 * t125;
t240 = pkin(4) * t214 + qJ(5) * t218;
t53 = pkin(4) * t203 + qJ(5) * t205;
t238 = t120 * t46;
t237 = t120 * t79;
t236 = t121 * t36;
t235 = t121 * t48;
t234 = t121 * t50;
t233 = t123 * t14;
t232 = t123 * t48;
t229 = t15 * t120;
t100 = t125 * qJ(2);
t227 = t122 * t98 + t100;
t225 = t125 * pkin(1) + t122 * qJ(2);
t224 = qJD(4) * t48;
t223 = qJD(4) * t50;
t128 = qJD(1) ^ 2;
t220 = t115 * t128;
t217 = t121 * t128;
t215 = t123 * t124;
t213 = t25 * qJD(4);
t212 = t26 * qJD(1);
t209 = t108 - t107;
t198 = qJD(4) * t124;
t196 = qJD(6) * t123;
t193 = qJDD(1) * t119;
t190 = qJDD(4) * t121;
t187 = g(1) * t214 + g(2) * t216 - t105;
t186 = t125 * qJ(3) + t225;
t185 = t124 * t237;
t184 = t79 * t215;
t183 = pkin(4) * t198 + qJ(5) * t200 + qJD(3);
t182 = t79 * t201;
t181 = t79 * t199;
t176 = qJDD(2) - t209;
t173 = qJD(4) * t242;
t172 = t207 * t51;
t171 = pkin(4) * t218 + t186;
t168 = qJD(1) * t207;
t127 = qJD(4) ^ 2;
t62 = qJDD(4) * t124 - t127 * t121;
t60 = t207 * qJDD(1);
t165 = 0.2e1 * t113 + t101 - t208;
t164 = -t119 - t103;
t163 = -t226 - t247;
t162 = t121 * t177;
t161 = -t117 + t176;
t149 = -t119 * t122 + t100;
t147 = t118 * t127 + t107;
t35 = -t163 + t119;
t59 = t242 * t124;
t19 = t120 * t59 + t123 * t35;
t18 = -t120 * t35 + t123 * t59;
t144 = t69 * qJD(3) - t119 * t160;
t141 = -t110 + t161;
t139 = t79 * t196 - t238;
t137 = -t111 + t141;
t133 = t26 * t203 + t166 + t187;
t11 = (-qJD(1) * qJD(5) - t195) * t124 + t134;
t29 = -qJD(5) * t124 + t183;
t132 = -qJD(1) * t29 - qJDD(1) * t57 - t11 + t147;
t131 = -t147 + t193 - t160 + t111;
t10 = -pkin(5) * t189 + (qJD(5) - t31) * qJD(4) - t167;
t129 = -qJD(6) * t126 * t79 - t121 * t208 + t10 - t246;
t95 = t116 * t128;
t84 = g(1) * t219;
t74 = t124 * t217;
t64 = t95 - t220;
t63 = t95 + t220;
t61 = t124 * t127 + t190;
t58 = t242 * t121;
t55 = 0.2e1 * t177 + t189;
t54 = -t94 + 0.2e1 * t178;
t45 = qJDD(1) * t116 - 0.2e1 * t162;
t44 = qJDD(1) * t115 + 0.2e1 * t162;
t43 = t62 - t217;
t42 = t190 + (t127 + t128) * t124;
t41 = t120 * t216 - t123 * t125;
t40 = t120 * t125 + t122 * t215;
t39 = t120 * t214 + t122 * t123;
t38 = t120 * t122 - t123 * t214;
t32 = pkin(8) * t203 + t53;
t28 = t121 * qJD(2) - t124 * t173;
t27 = -t124 * qJD(2) - t121 * t173;
t24 = -0.2e1 * t121 * t94 + 0.2e1 * (t115 - t116) * t194;
t23 = t124 * t170 + t183;
t13 = t120 * t30 + t123 * t32;
t12 = -t120 * t32 + t123 * t30;
t4 = -qJD(6) * t19 - t120 * t23 + t123 * t27;
t3 = qJD(6) * t18 + t120 * t27 + t123 * t23;
t8 = [0, 0, 0, 0, 0, qJDD(1), t209, t208, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t117 + t176, t165, t188 * pkin(1) - g(1) * (-t122 * pkin(1) + t100) - g(2) * t225 + (t101 + t113) * qJ(2), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) + t165, 0.2e1 * t111 - t141 + t193, -g(1) * t149 - g(2) * t186 + qJ(2) * t175 + t92 * qJD(2) + t144, t45, t24, t62, t44, -t61, 0, t131 * t121 + t124 * t254 + t84, -t254 * t121 + (t131 + t108) * t124, t208 + t207 * (-qJDD(1) * t118 - t112 - t51) -g(1) * (t149 - t245) - g(2) * (-t122 * pkin(7) + t186) + t118 * t172 + t207 * t68 * qJD(2) + t144, 0, -t62, t61, t45, t24, t44, -qJD(2) * t168 - t118 * t60 - t261, t132 * t121 - t124 * t253 - t84, t253 * t121 + (t132 - t108) * t124, t11 * t57 + t26 * t29 - g(1) * (t227 - t245) - g(2) * (-qJ(5) * t214 + t171) + (g(2) * pkin(7) - g(1) * t164) * t122 + (-t124 * t33 - t236) * qJD(2) + t130 * t118, t196 * t234 + (-t121 * t14 + t198 * t50) * t120 (-t120 * t48 + t123 * t50) * t198 + (-t229 - t233 + (-t120 * t50 - t232) * qJD(6)) * t121 (-t14 + t182) * t124 + (t139 - t223) * t121, -t198 * t232 + (-t123 * t15 + t197 * t48) * t121 (-t15 + t181) * t124 + (t138 + t224) * t121, -t124 * t46 - t200 * t79, -g(1) * t41 + g(2) * t39 - t15 * t58 - t18 * t46 + t28 * t48 + t4 * t79 + (-t199 * t25 + t2) * t124 + (-qJD(4) * t5 - t10 * t123 + t197 * t25) * t121, -g(1) * t40 - g(2) * t38 + t14 * t58 + t19 * t46 + t28 * t50 - t3 * t79 + (t201 * t25 - t1) * t124 + (qJD(4) * t6 + t10 * t120 + t196 * t25) * t121, t14 * t18 - t15 * t19 - t3 * t48 - t4 * t50 + t84 - t155 * t198 + (-t107 - t251) * t121, t1 * t19 + t6 * t3 + t2 * t18 + t5 * t4 - t10 * t58 + t25 * t28 - g(1) * t227 - g(2) * t171 + (g(1) * t250 - g(2) * t148) * t125 + (-g(1) * (t164 - t247) + g(2) * t250) * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t128, -qJ(2) * t128 + t161, 0, 0, 0, 0, 0, 0, 0, -t128, -qJDD(1), -qJD(1) * t92 + t137, 0, 0, 0, 0, 0, 0, -t55, t54, t63, -t168 * t68 + t137, 0, 0, 0, 0, 0, 0, t63, t55, -t54, qJ(5) * t94 + (t236 + (qJD(5) + t33) * t124) * qJD(1) - t134 - t209, 0, 0, 0, 0, 0, 0 (t184 - t235) * qJD(1) + t139 (-t185 - t234) * qJD(1) + t138, t135 * t120 + t262 * t123 (-t121 * t25 + t124 * t156) * qJD(1) - t209 + t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t128, -t263 + t175, 0, 0, 0, 0, 0, 0, t43, -t42, -t60, t172 - t263, 0, 0, 0, 0, 0, 0, -t60, -t43, t42, -t212 + t261, 0, 0, 0, 0, 0, 0, qJD(1) * t237 + (t15 + t181) * t121 + (-t138 + t224) * t124, t79 * t204 + (-t14 - t182) * t121 + (t139 + t223) * t124 (-t50 * t200 + qJD(1) * t48 + (qJD(6) * t48 - t14) * t124) * t123 + (-t48 * t200 - qJD(1) * t50 + (t15 - t222) * t124) * t120, t155 * qJD(1) + (qJD(4) * t156 + t10) * t121 + (t213 - t252) * t124 - t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t64, t94, -t74, -t189, qJDD(4) (t51 - t211) * t124 - t187, t263 * t121 + t246 - t37, 0, 0, qJDD(4), -t94, t189, t74, t64, -t74 (-pkin(4) * t124 - qJ(5) * t121) * qJDD(1) + ((-t36 - t202) * t124 + (-qJD(5) + t33 + t239) * t121) * qJD(1), t205 * t53 + t133 - 0.2e1 * t221, 0.2e1 * t192 + 0.2e1 * qJD(4) * qJD(5) + t37 + (qJD(1) * t53 - g(3)) * t124 + (-t208 - t212) * t121, -t17 * pkin(4) - g(1) * t240 - g(2) * t241 + g(3) * t226 - t16 * qJ(5) + t169 * t36 - t228 * t68 - t26 * t53, -t237 * t50 - t233 (-t15 - t243) * t123 + (t14 + t142) * t120 (-t185 + t234) * qJD(1) + t138, t123 * t142 + t229 (-t184 - t235) * qJD(1) - t139, t79 * t205, qJ(5) * t15 - t12 * t79 + t129 * t120 + t123 * t255 + t5 * t205 + t210 * t48, -qJ(5) * t14 - t120 * t255 + t129 * t123 + t13 * t79 - t6 * t205 + t210 * t50, t12 * t50 + t13 * t48 + (-t6 * t203 + t126 * t14 - t2 + (-t126 * t48 - t6) * qJD(6)) * t123 + (t5 * t203 - t126 * t15 - t1 + (t126 * t50 + t5) * qJD(6)) * t120 - t187, t10 * qJ(5) - t6 * t13 - t5 * t12 - g(1) * (pkin(8) * t214 + t240) - g(2) * (pkin(8) * t216 + t241) - g(3) * t163 + t210 * t25 + t252 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, qJDD(4) - t74, -t95 - t127, qJD(4) * t36 + t133 - t221 + t47, 0, 0, 0, 0, 0, 0, -t237 * t79 - t224 - t34, -t123 * t79 ^ 2 - t223 + t238, -t120 * t262 + t135 * t123, -t213 + t259 * t123 + (t1 - t249) * t120 + t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, -t48 ^ 2 + t50 ^ 2, -t135, -t244, -t262, -t46, -g(1) * t38 + g(2) * t40 - t105 * t123 - t25 * t50 + t259, -g(1) * t39 - g(2) * t41 + t25 * t48 + t249 + (-qJD(6) * t21 - t7) * t123 + (qJD(6) * t20 + t105 - t9) * t120, 0, 0;];
tau_reg  = t8;
