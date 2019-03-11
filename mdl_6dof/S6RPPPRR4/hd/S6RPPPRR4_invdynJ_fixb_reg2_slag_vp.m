% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPPRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:59
% EndTime: 2019-03-09 01:36:02
% DurationCPUTime: 2.50s
% Computational Cost: add. (3346->377), mult. (5381->456), div. (0->0), fcn. (3075->8), ass. (0->207)
t107 = cos(qJ(5));
t178 = t107 * qJDD(1);
t253 = qJD(5) * qJD(6) - t178;
t105 = sin(qJ(5));
t101 = sin(pkin(9));
t102 = cos(pkin(9));
t200 = qJ(2) * qJD(1);
t108 = -pkin(1) - pkin(2);
t71 = t108 * qJD(1) + qJD(2);
t45 = -t101 * t200 + t102 * t71;
t40 = qJD(1) * pkin(3) + qJD(4) - t45;
t37 = qJD(1) * pkin(7) + t40;
t29 = qJD(3) * t107 + t105 * t37;
t208 = qJD(5) * t29;
t249 = qJDD(1) * pkin(3) + qJDD(4);
t183 = qJDD(1) * t101;
t70 = t108 * qJDD(1) + qJDD(2);
t226 = qJ(2) * t183 - t102 * t70;
t187 = qJD(1) * qJD(2);
t76 = t101 * t187;
t38 = -t76 - t226;
t34 = -t38 + t249;
t30 = qJDD(1) * pkin(7) + t34;
t23 = t107 * t30;
t10 = -t105 * qJDD(3) - t208 + t23;
t195 = qJD(5) * t107;
t177 = -t107 * qJDD(3) - t105 * t30 - t37 * t195;
t185 = qJD(3) * qJD(5);
t9 = -t105 * t185 - t177;
t145 = t10 * t107 + t9 * t105;
t252 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4);
t237 = sin(qJ(1));
t238 = cos(qJ(1));
t54 = -t237 * t101 - t238 * t102;
t241 = g(2) * t54;
t55 = t238 * t101 - t237 * t102;
t242 = g(1) * t55;
t154 = -t241 + t242;
t8 = -qJDD(5) * pkin(5) - t10;
t119 = -g(3) * t105 - t107 * t154 - t8;
t74 = qJD(1) * t105 - qJD(6);
t251 = pkin(8) * qJD(6) * t74 + t119;
t110 = qJD(1) ^ 2;
t182 = qJDD(1) * t102;
t66 = t101 * t110 + t182;
t197 = qJD(5) * t105;
t28 = -qJD(3) * t105 + t107 * t37;
t20 = -qJD(5) * pkin(5) - t28;
t250 = -t107 * t8 + t20 * t197;
t104 = sin(qJ(6));
t106 = cos(qJ(6));
t21 = qJD(5) * pkin(8) + t29;
t148 = pkin(5) * t105 - pkin(8) * t107;
t211 = qJ(2) * t102;
t128 = -t148 + t211;
t221 = t101 * t71;
t96 = qJD(1) * qJ(4);
t31 = qJD(1) * t128 + t221 - t96;
t11 = -t104 * t21 + t106 * t31;
t12 = t104 * t31 + t106 * t21;
t142 = t104 * t11 - t106 * t12;
t248 = qJD(5) * t142 + t8;
t180 = qJDD(5) * t107;
t109 = qJD(5) ^ 2;
t201 = t109 + t110;
t247 = t201 * t105 - t180;
t181 = qJDD(5) * t105;
t246 = t201 * t107 + t181;
t198 = qJD(2) * t101;
t46 = t102 * t200 + t221;
t41 = t46 - t96;
t63 = -t101 * qJ(2) + t102 * t108;
t59 = pkin(3) - t63;
t52 = pkin(7) + t59;
t205 = t101 * t108;
t64 = t205 + t211;
t56 = -qJ(4) + t64;
t245 = qJD(5) * (qJD(1) * t56 - t198 + t41) - qJDD(5) * t52;
t244 = -(t105 * t28 - t107 * t29) * qJD(5) + t145;
t243 = pkin(7) * t54;
t240 = g(2) * t55;
t239 = t55 * pkin(7);
t236 = g(3) * t107;
t234 = t11 * t74;
t233 = t12 * t74;
t196 = qJD(5) * t106;
t199 = qJD(1) * t107;
t57 = t104 * t199 + t196;
t189 = t104 * qJD(5);
t58 = t106 * t199 - t189;
t232 = t57 * t58;
t231 = t57 * t74;
t230 = t58 * t74;
t167 = t102 * t195;
t202 = t105 * t106;
t203 = t104 * t105;
t48 = t101 * t106 + t102 * t203;
t229 = qJD(6) * t48 - t106 * t167 - (t101 * t202 + t102 * t104) * qJD(1);
t132 = -t101 * t104 + t102 * t202;
t228 = qJD(6) * t132 + t104 * t167 - (-t101 * t203 + t102 * t106) * qJD(1);
t227 = qJ(2) * t182 + t101 * t70;
t225 = t238 * pkin(1) + t237 * qJ(2);
t224 = g(1) * t237 - g(2) * t238;
t97 = t105 ^ 2;
t98 = t107 ^ 2;
t223 = t97 - t98;
t222 = t97 + t98;
t220 = t104 * t54;
t219 = t104 * t57;
t218 = t104 * t58;
t217 = t106 * t54;
t216 = t106 * t57;
t215 = t106 * t58;
t192 = qJD(6) * t107;
t166 = t104 * t192;
t24 = t104 * qJDD(5) + (t105 * t196 + t166) * qJD(1) + t253 * t106;
t214 = t24 * t104;
t186 = qJD(1) * qJD(5);
t164 = t105 * t186;
t25 = (qJD(1) * t192 + qJDD(5)) * t106 + (-t164 - t253) * t104;
t213 = t25 * t106;
t212 = pkin(1) * qJDD(1);
t210 = qJD(1) * t41;
t209 = qJD(1) * t74;
t207 = qJD(5) * t57;
t206 = qJD(5) * t58;
t99 = qJDD(3) + g(3);
t194 = qJD(6) * t104;
t193 = qJD(6) * t106;
t190 = t102 * qJD(2);
t188 = qJ(2) * qJDD(1);
t179 = t105 * qJDD(1);
t176 = t238 * pkin(2) + t225;
t174 = 0.2e1 * t187;
t173 = t74 * t189;
t171 = t74 * t196;
t170 = t107 * t110 * t105;
t168 = t101 * t199;
t165 = t102 * t187;
t163 = t107 * t186;
t162 = -t54 * pkin(3) + t176;
t161 = -t227 - t252;
t160 = -qJD(6) * t57 + t24;
t159 = -qJD(6) * t58 + t25;
t157 = qJDD(2) - t212;
t156 = t105 * t163;
t155 = -t237 * pkin(1) + t238 * qJ(2);
t153 = g(1) * t54 + t240;
t152 = t25 + t173;
t151 = -t24 + t171;
t149 = -pkin(5) * t107 - pkin(8) * t105;
t118 = qJD(5) * t149 + t190;
t150 = -qJD(6) * t105 * t52 - qJD(4) + t118;
t39 = t165 + t227;
t32 = t39 + t252;
t73 = -qJD(4) + t190;
t147 = t32 * t56 + t41 * t73;
t146 = -qJD(6) * t21 + t240;
t144 = t101 * t45 - t102 * t46;
t143 = t104 * t12 + t106 * t11;
t141 = t105 * t29 + t107 * t28;
t139 = qJ(4) + t148;
t138 = qJ(4) * t55 + t162;
t137 = t107 * t24 + t58 * t197;
t136 = -t107 * t25 + t57 * t197;
t53 = -qJDD(6) + t163 + t179;
t135 = -t104 * t53 - t74 * t193;
t134 = -t106 * t53 + t74 * t194;
t133 = -t154 + t226;
t130 = g(1) * t238 + g(2) * t237;
t129 = -qJDD(1) * t56 - t153;
t127 = -t237 * pkin(2) + t155;
t126 = pkin(8) * t53 - t20 * t74;
t125 = t133 + t249;
t124 = t55 * pkin(3) + t127;
t123 = -t135 - t207;
t122 = -t134 - t206;
t121 = -t105 * t154 + t236;
t7 = qJDD(5) * pkin(8) + t9;
t120 = -qJD(6) * t31 - t105 * t241 - t236 - t7;
t117 = t154 + t185 - t210;
t42 = -qJ(4) + t128 + t205;
t116 = qJD(6) * t42 + t105 * t198 + t52 * t195;
t115 = t29 * t195 - t28 * t197 + t145 - t154;
t114 = t54 * qJ(4) + t124;
t14 = qJD(1) * t118 - qJDD(1) * t148 - t161;
t1 = qJD(6) * t11 + t104 * t14 + t106 * t7;
t13 = t106 * t14;
t2 = -qJD(6) * t12 - t104 * t7 + t13;
t113 = -qJD(6) * t143 + t1 * t106 - t2 * t104;
t112 = -qJD(1) * t73 - t109 * t52 + t129 - t32;
t111 = qJD(5) * t20 + t113;
t85 = t98 * qJDD(1);
t84 = t97 * qJDD(1);
t68 = t105 * t109 - t180;
t67 = t107 * t109 + t181;
t65 = -t102 * t110 + t183;
t62 = t149 * qJD(1);
t27 = t55 * t202 - t220;
t26 = -t55 * t203 - t217;
t18 = t104 * t42 + t52 * t202;
t17 = t106 * t42 - t52 * t203;
t16 = t104 * t62 + t106 * t28;
t15 = -t104 * t28 + t106 * t62;
t4 = -t104 * t116 + t106 * t150;
t3 = t104 * t150 + t106 * t116;
t5 = [0, 0, 0, 0, 0, qJDD(1), t224, t130, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t212 + t224, 0, -t130 + t174 + 0.2e1 * t188, -t157 * pkin(1) - g(1) * t155 - g(2) * t225 + (t174 + t188) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), -qJDD(1) * t63 + t133 + 0.2e1 * t76, qJDD(1) * t64 + t153 + 0.2e1 * t165 + t227, 0, -g(1) * t127 - g(2) * t176 - t144 * qJD(2) + t38 * t63 + t39 * t64, qJDD(1), 0, 0, 0, 0, 0, 0, -qJDD(1) * t59 - t125 - 0.2e1 * t76 (-t73 - t190) * qJD(1) + t129 + t161, -g(1) * t114 - g(2) * t138 + t40 * t198 + t34 * t59 + t147, t85 - 0.2e1 * t156, -0.2e1 * t105 * t178 + 0.2e1 * t223 * t186, t68, t84 + 0.2e1 * t156, t67, 0, t112 * t105 - t245 * t107, t245 * t105 + t112 * t107, t115 + (qJDD(1) * t52 + t76) * t222, -g(1) * (t114 + t239) - g(2) * (t138 - t243) + t141 * t198 + t244 * t52 + t147, -t106 * t137 - t166 * t58 (t216 + t218) * t197 + (t214 - t213 + (-t215 + t219) * qJD(6)) * t107 (-t24 - t171) * t105 + (-t134 + t206) * t107, -t104 * t136 + t192 * t216 (-t25 + t173) * t105 + (t135 - t207) * t107, t105 * t53 + t195 * t74, -t104 * t242 - g(2) * t27 - t17 * t53 - t4 * t74 + (-g(1) * t217 - t2 + (t104 * t20 - t52 * t57) * qJD(5)) * t105 + (-qJD(5) * t11 - t104 * t8 - t193 * t20 + t198 * t57 + t25 * t52) * t107, -t106 * t242 - g(2) * t26 + t18 * t53 + t3 * t74 + (g(1) * t220 + t1 + (t106 * t20 - t52 * t58) * qJD(5)) * t105 + (qJD(5) * t12 - t106 * t8 + t194 * t20 + t198 * t58 - t24 * t52) * t107, -t17 * t24 + t18 * t25 + t3 * t57 + t4 * t58 - t143 * t197 + (-qJD(6) * t142 + t1 * t104 + t106 * t2 + t153) * t107, t1 * t18 + t12 * t3 + t2 * t17 + t11 * t4 - t20 * t107 * t198 - g(1) * (t139 * t54 + t124 + t239) - g(2) * (t139 * t55 + t162 - t243) + t250 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t110, -qJ(2) * t110 + t157 - t224, 0, 0, 0, 0, 0, 0, -t66, t65, 0, qJD(1) * t144 + t101 * t39 + t102 * t38 - t224, 0, 0, 0, 0, 0, 0, 0, t66, -t65, t101 * t32 - t102 * t34 + (-t101 * t40 - t102 * t41) * qJD(1) - t224, 0, 0, 0, 0, 0, 0 (-0.2e1 * t163 - t179) * t101 + t247 * t102 (0.2e1 * t164 - t178) * t101 + t246 * t102, -t66 * t222 (-qJD(1) * t141 + t32) * t101 + (-t210 - t244) * t102 - t224, 0, 0, 0, 0, 0, 0, t102 * t136 - t168 * t57 - t228 * t74 - t48 * t53, t102 * t137 - t132 * t53 - t168 * t58 + t229 * t74, -t132 * t25 + t228 * t58 + t229 * t57 - t24 * t48, -t1 * t132 - t250 * t102 + t228 * t11 + t229 * t12 + t20 * t168 + t2 * t48 - t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, 0, 0, 0, 0, 0, -t67, t68, 0, -qJD(5) * t141 - t10 * t105 + t9 * t107 + g(3), 0, 0, 0, 0, 0, 0, -t105 * t152 + t107 * t123, -t105 * t151 + t107 * t122 (-t216 + t218) * t197 + (t214 + t213 + (-t215 - t219) * qJD(6)) * t107, t248 * t105 + t111 * t107 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), -t110, t125 + t76 + t210, 0, 0, 0, 0, 0, 0, -t247, -t246, t84 + t85, t115 + t210, 0, 0, 0, 0, 0, 0, t105 * t123 - t106 * t209 + t107 * t152, t104 * t209 + t105 * t122 + t107 * t151 (qJD(1) * t58 + t105 * t159 + t195 * t57) * t106 + (qJD(1) * t57 + t105 * t160 - t195 * t58) * t104, t143 * qJD(1) + t111 * t105 - t248 * t107 - t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, -t223 * t110, -t178, -t170, t179, qJDD(5), t208 + t23 + (-qJD(5) * t37 - t99) * t105 - t117 * t107, qJD(5) * t28 + t105 * t117 + t177 - t236, 0, 0, t215 * t74 + t214 (t24 - t231) * t106 + (t25 - t230) * t104 (-t107 * t58 + t202 * t74) * qJD(1) + t135, t219 * t74 + t213 (t107 * t57 - t203 * t74) * qJD(1) + t134, -t74 * t199, pkin(5) * t25 + t126 * t104 + t251 * t106 + t11 * t199 + t15 * t74 + t29 * t57, -pkin(5) * t24 - t251 * t104 + t126 * t106 - t12 * t199 - t16 * t74 + t29 * t58, -t15 * t58 - t16 * t57 + (pkin(8) * t159 + t1 + t234) * t106 + (pkin(8) * t160 - t2 + t233) * t104 + t121, -t11 * t15 - t12 * t16 - t20 * t29 + t119 * pkin(5) + (t113 + t121) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, -t57 ^ 2 + t58 ^ 2, t24 + t231, -t232, t25 + t230, -t53, -g(1) * t26 + t104 * t120 + t106 * t146 + t20 * t58 + t13 - t233, g(1) * t27 - t234 - t20 * t57 + (-t14 - t146) * t104 + t120 * t106, 0, 0;];
tau_reg  = t5;
