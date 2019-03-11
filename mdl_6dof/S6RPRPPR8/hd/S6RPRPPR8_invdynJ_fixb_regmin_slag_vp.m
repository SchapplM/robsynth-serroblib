% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:00:03
% EndTime: 2019-03-09 03:00:08
% DurationCPUTime: 2.38s
% Computational Cost: add. (1494->346), mult. (2623->398), div. (0->0), fcn. (1383->6), ass. (0->202)
t110 = cos(qJ(3));
t192 = qJD(1) * t110;
t184 = qJD(1) * qJD(3);
t165 = t110 * t184;
t107 = sin(qJ(3));
t178 = t107 * qJDD(1);
t246 = t165 + t178;
t112 = -pkin(3) - pkin(4);
t113 = -pkin(1) - pkin(7);
t59 = t113 * qJD(1) + qJD(2);
t30 = (qJ(5) * qJD(1) + t59) * t110;
t245 = qJD(4) - t30;
t20 = qJD(3) * t112 + t245;
t100 = qJD(3) * qJD(4);
t189 = qJD(3) * t110;
t47 = t59 * t189;
t240 = t113 * qJDD(1);
t58 = qJDD(2) + t240;
t49 = t107 * t58;
t99 = qJDD(3) * qJ(4);
t15 = t47 + t49 + t99 + t100;
t191 = qJD(3) * t107;
t46 = t59 * t191;
t50 = t110 * t58;
t170 = qJDD(4) + t46 - t50;
t209 = qJDD(3) * pkin(3);
t16 = t170 - t209;
t155 = -t110 * t59 + qJD(4);
t224 = qJD(3) * pkin(3);
t32 = t155 - t224;
t102 = qJD(3) * qJ(4);
t51 = t107 * t59;
t34 = t51 + t102;
t117 = (t107 * t32 + t110 * t34) * qJD(3) + t15 * t107 - t16 * t110;
t111 = cos(qJ(1));
t96 = g(2) * t111;
t108 = sin(qJ(1));
t97 = g(1) * t108;
t243 = t97 - t96;
t244 = t117 - t243;
t242 = g(1) * t111 + g(2) * t108;
t201 = qJ(5) + t113;
t167 = qJDD(3) * t112;
t183 = qJD(1) * qJD(5);
t239 = -t246 * qJ(5) - t107 * t183;
t204 = t110 * t111;
t206 = t108 * t110;
t94 = g(3) * t107;
t238 = g(1) * t206 - g(2) * t204 - t94;
t101 = -pkin(8) + t112;
t185 = qJ(5) * qJDD(1);
t166 = t107 * t184;
t62 = qJ(5) * t166;
t118 = t62 + (-t183 - t185) * t110 + t170;
t129 = t110 * pkin(5) + t101 * t107 - qJ(2);
t82 = qJ(4) * t192;
t198 = qJD(5) + t82;
t14 = qJD(1) * t129 + t198;
t162 = -qJD(6) * t14 - t101 * qJDD(3) - t118;
t194 = qJD(1) * t107;
t81 = qJ(5) * t194;
t25 = -t81 - t34;
t21 = qJD(3) * pkin(5) - t25;
t88 = t110 * qJ(4);
t23 = t88 + t129;
t153 = qJD(3) * t201;
t26 = -t110 * qJD(5) + t107 * t153;
t84 = t110 * qJDD(1);
t42 = -qJDD(6) - t84 + t166;
t53 = t201 * t110;
t9 = -t15 + t239;
t6 = qJDD(3) * pkin(5) - t9;
t67 = qJD(6) + t192;
t237 = -(qJD(6) * t23 + t26) * t67 - t53 * t42 + t6 * t107 + (qJD(3) * t21 + t162) * t110;
t213 = qJ(4) * t107;
t133 = t110 * t112 - t213;
t106 = sin(qJ(6));
t109 = cos(qJ(6));
t186 = t106 * qJD(3);
t187 = qJD(6) * t109;
t12 = (t107 * t187 + t110 * t186) * qJD(1) - qJD(6) * t186 + t109 * qJDD(3) + t106 * t178;
t105 = qJ(4) + pkin(5);
t127 = t101 * t110 - t105 * t107;
t234 = g(3) * t110;
t236 = t107 * t243 + (t127 * qJD(1) + qJD(6) * t101) * t67 - t6 + t234;
t235 = pkin(3) * t107;
t233 = t23 * t42;
t190 = qJD(3) * t109;
t43 = t106 * t194 + t190;
t232 = t43 * t67;
t193 = qJD(1) * t109;
t45 = t107 * t193 - t186;
t231 = t45 * t67;
t208 = t107 * t108;
t230 = pkin(3) * t206 + qJ(4) * t208;
t229 = g(1) * t204 + g(2) * t206;
t86 = t110 * qJD(4);
t228 = -qJ(4) * t84 - qJD(1) * t86;
t227 = t111 * pkin(1) + t108 * qJ(2);
t223 = t106 * t42;
t222 = t106 * t67;
t221 = t107 * t21;
t220 = t107 * t43;
t219 = t107 * t45;
t218 = t109 * t42;
t217 = t109 * t45;
t11 = -qJD(6) * t43 - t106 * qJDD(3) + t246 * t109;
t216 = t11 * t106;
t215 = pkin(1) * qJDD(1);
t115 = qJD(1) ^ 2;
t214 = qJ(2) * t115;
t212 = qJD(3) * t25;
t211 = qJD(3) * t43;
t210 = qJD(3) * t45;
t207 = t107 * t111;
t205 = t109 * t110;
t203 = t110 * t115;
t114 = qJD(3) ^ 2;
t202 = t113 * t114;
t173 = t112 * t107;
t148 = -qJ(2) + t173;
t24 = qJD(1) * t148 + t198;
t199 = -qJD(5) - t24;
t103 = t107 ^ 2;
t104 = t110 ^ 2;
t197 = -t103 - t104;
t196 = t103 - t104;
t195 = t114 + t115;
t188 = qJD(6) * t106;
t182 = qJDD(1) * qJ(2);
t181 = qJDD(3) * t107;
t180 = qJDD(3) * t110;
t179 = qJDD(3) * t113;
t177 = t20 * t191 - t243;
t176 = t110 * t222;
t175 = t67 * t205;
t174 = 0.2e1 * qJD(1) * qJD(2);
t172 = t67 * t186;
t171 = t67 * t190;
t169 = qJDD(5) - t228;
t164 = qJ(2) + t235;
t8 = t167 + t118;
t163 = -t8 - t212;
t19 = t101 * qJD(3) + t245;
t119 = qJD(3) * t127 - qJD(2);
t2 = qJD(1) * t119 + qJDD(1) * t129 + t169;
t161 = qJD(6) * t19 - t2;
t160 = pkin(3) * t208 + t111 * pkin(7) + t227;
t159 = -t50 + t238;
t41 = t88 + t148;
t158 = qJD(1) * t41 + t24;
t152 = -0.2e1 * t166;
t151 = 0.2e1 * t165;
t150 = qJDD(2) - t215;
t149 = -g(2) * t207 + t234 - t49;
t145 = -qJDD(4) - t159;
t144 = pkin(3) * t110 + t213;
t33 = t164 * qJD(1) - t82;
t54 = t164 - t88;
t141 = (qJD(1) * t54 + t33) * qJD(3);
t89 = t111 * qJ(2);
t140 = pkin(3) * t207 - qJ(4) * t204 + t89;
t139 = t109 * t67;
t138 = t162 + t94;
t137 = t174 + 0.2e1 * t182;
t136 = -t67 * t187 + t223;
t135 = t67 * t188 + t218;
t122 = qJD(3) * t133 - qJD(2);
t22 = t86 + t122;
t121 = qJDD(1) * t148 + t169;
t7 = qJD(1) * t122 + t121;
t134 = qJD(1) * t22 + qJDD(1) * t41 + t7;
t132 = t46 + t62 - t145;
t126 = qJD(3) * t144 + qJD(2);
t125 = t137 - t202;
t10 = t126 * qJD(1) + t164 * qJDD(1) + t228;
t28 = -t86 + t126;
t124 = -qJD(1) * t28 - qJDD(1) * t54 - t10 + t202;
t123 = -g(1) * t208 + 0.2e1 * t100 - t149 + 0.2e1 * t99;
t29 = t51 + t81;
t120 = t101 * t42 + (-t21 + t29) * t67;
t83 = t104 * qJDD(1);
t65 = t107 * t203;
t64 = t110 * t179;
t60 = -t104 * t115 - t114;
t56 = qJDD(3) - t65;
t55 = -qJDD(1) * t103 - t83;
t52 = t201 * t107;
t48 = t144 * qJD(1);
t40 = -t195 * t107 + t180;
t39 = t195 * t110 + t181;
t38 = t106 * t108 - t109 * t204;
t37 = t106 * t204 + t108 * t109;
t36 = t106 * t111 + t108 * t205;
t35 = t106 * t206 - t109 * t111;
t31 = t133 * qJD(1);
t27 = t107 * qJD(5) + t110 * t153;
t13 = t86 + t119;
t4 = t106 * t14 + t109 * t19;
t3 = -t106 * t19 + t109 * t14;
t1 = t109 * t2;
t5 = [qJDD(1), t243, t242, qJDD(2) - 0.2e1 * t215 - t243, t137 - t242, -t150 * pkin(1) - g(1) * (-t108 * pkin(1) + t89) - g(2) * t227 + (t174 + t182) * qJ(2), t110 * t152 + t83, -0.2e1 * t107 * t84 + 0.2e1 * t196 * t184, -t107 * t114 + t180, -t110 * t114 - t181, 0, qJ(2) * t151 + t64 + (t125 - t242) * t107 (-0.2e1 * qJ(2) * t184 - t179) * t107 + t125 * t110 - t229, t64 + t110 * t141 + (-t124 - t242) * t107, t197 * t240 - t244 (t141 + t179) * t107 + t124 * t110 + t229, t10 * t54 + t33 * t28 - g(1) * (t113 * t108 + t140) - g(2) * (-t108 * t88 + t160) + t117 * t113, t52 * qJDD(3) + t134 * t110 + (-t107 * t158 + t27) * qJD(3) + t229, -t53 * qJDD(3) + (t110 * t158 + t26) * qJD(3) + (t134 + t242) * t107 (qJDD(1) * t52 - t9 + (-qJD(3) * t53 + t27) * qJD(1)) * t107 + (qJDD(1) * t53 + (qJD(3) * t52 - t26) * qJD(1) + t163) * t110 + t177, -t8 * t53 + t20 * t26 - t9 * t52 - t25 * t27 + t7 * t41 + t24 * t22 - g(1) * (pkin(4) * t207 + t140) - g(2) * (-t111 * qJ(5) + t160) + (-g(1) * t201 - g(2) * (pkin(4) * t107 - t88)) * t108, t189 * t217 + (t11 * t109 - t45 * t188) * t107 (-t106 * t45 - t109 * t43) * t189 + (-t216 - t109 * t12 + (t106 * t43 - t217) * qJD(6)) * t107 (t11 + t171) * t110 + (-t135 - t210) * t107 (-t12 - t172) * t110 + (t136 + t211) * t107, -t110 * t42 - t67 * t191, -t3 * t191 - g(1) * t38 + g(2) * t36 + t1 * t110 + t52 * t12 + t27 * t43 + (t13 * t67 - t233 + (-t110 * t19 + t53 * t67 + t221) * qJD(6)) * t109 + t237 * t106, t4 * t191 - g(1) * t37 - g(2) * t35 + t52 * t11 + t27 * t45 + (-(qJD(6) * t53 + t13) * t67 + t233 + t161 * t110 - qJD(6) * t221) * t106 + t237 * t109; 0, 0, 0, qJDD(1), -t115, t150 - t214 - t243, 0, 0, 0, 0, 0, t40, -t39, t40, t55, t39, -t33 * qJD(1) + t244, t39, -t40, -t55, t24 * qJD(1) - t9 * t107 + t163 * t110 + t177, 0, 0, 0, 0, 0, t67 * t193 + (t12 - t172) * t107 + (-t136 + t211) * t110, -qJD(1) * t222 + (t11 - t171) * t107 + (-t135 + t210) * t110; 0, 0, 0, 0, 0, 0, t65, -t196 * t115, t84, -t178, qJDD(3), -qJ(2) * t203 - t159 (t214 + t97) * t107 + t149, 0.2e1 * t209 + (-t107 * t48 - t110 * t33) * qJD(1) + t145, -t144 * qJDD(1) + ((t34 - t102) * t110 + (-qJD(4) + t32 + t224) * t107) * qJD(1) (-t107 * t33 + t110 * t48) * qJD(1) + t123, t15 * qJ(4) - t16 * pkin(3) - t33 * t48 - t32 * t51 - g(1) * t230 - g(3) * (t88 - t235) + t155 * t34 + t144 * t96, -qJD(3) * t30 + t47 + (t107 * t24 - t110 * t31) * qJD(1) + t123 - t239, -qJ(5) * t84 - qJD(3) * t29 + 0.2e1 * t167 + (-t107 * t31 + t199 * t110) * qJD(1) + t132, -t133 * qJDD(1) + (t25 + t29 + t102) * t192, t8 * t112 - t9 * qJ(4) - t20 * t29 - t24 * t31 - g(1) * (pkin(4) * t206 + t230) - g(3) * (t88 + t173) - t245 * t25 - t133 * t96, -t139 * t45 - t216 (-t11 + t232) * t109 + (t12 + t231) * t106 (-t175 + t219) * qJD(1) + t136 (t176 - t220) * qJD(1) + t135, t67 * t194, t105 * t12 + t120 * t106 - t109 * t236 + t3 * t194 + t245 * t43, t105 * t11 + t106 * t236 + t120 * t109 - t4 * t194 + t245 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t84, t60, -qJD(3) * t34 + t33 * t192 + t16 + t238, t60, t56, -t84, t212 + t167 + (t199 * qJD(1) - t185) * t110 + t132, 0, 0, 0, 0, 0, -t139 * t67 - t211 + t223, t106 * t67 ^ 2 - t210 + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 + t152, t151 + t178, t197 * t115 (t107 * t25 + t110 * t20 + t122) * qJD(1) + t121 + t242, 0, 0, 0, 0, 0 (-t176 - t220) * qJD(1) - t135 (-t175 - t219) * qJD(1) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t43, -t43 ^ 2 + t45 ^ 2, t11 + t232, -t12 + t231, -t42, -g(1) * t35 + g(2) * t37 + t106 * t138 - t187 * t19 - t21 * t45 + t4 * t67 + t1, -g(1) * t36 - g(2) * t38 + t106 * t161 + t109 * t138 + t21 * t43 + t3 * t67;];
tau_reg  = t5;
