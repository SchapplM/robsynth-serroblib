% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% tau_reg [6x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPPRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:37
% EndTime: 2019-03-08 18:40:45
% DurationCPUTime: 2.64s
% Computational Cost: add. (3171->305), mult. (8474->488), div. (0->0), fcn. (9030->18), ass. (0->184)
t115 = cos(qJ(5));
t194 = qJD(4) * t115;
t90 = -qJD(6) + t194;
t237 = t90 + qJD(6);
t113 = sin(qJ(4));
t116 = cos(qJ(4));
t103 = sin(pkin(8));
t108 = cos(pkin(8));
t100 = sin(pkin(14));
t101 = sin(pkin(13));
t105 = sin(pkin(6));
t106 = cos(pkin(14));
t107 = cos(pkin(13));
t109 = cos(pkin(7));
t197 = t107 * t109;
t142 = (-t100 * t101 + t106 * t197) * t105;
t104 = sin(pkin(7));
t201 = t104 * t106;
t110 = cos(pkin(6));
t89 = qJD(1) * t110 + qJD(2);
t55 = qJD(1) * t142 + t89 * t201;
t219 = t108 * t55;
t199 = t105 * t104;
t181 = t107 * t199;
t71 = -qJD(1) * t181 + t109 * t89 + qJD(3);
t165 = t103 * t71 + t219;
t206 = t101 * t106;
t141 = (t100 * t197 + t206) * t105;
t207 = t100 * t104;
t56 = qJD(1) * t141 + t89 * t207;
t140 = t113 * t56 - t165 * t116;
t86 = t110 * qJDD(1) + qJDD(2);
t53 = qJDD(1) * t142 + t86 * t201;
t70 = -qJDD(1) * t181 + t109 * t86 + qJDD(3);
t166 = t103 * t70 + t108 * t53;
t195 = qJD(4) * t113;
t179 = t103 * t195;
t193 = qJD(4) * t116;
t54 = qJDD(1) * t141 + t86 * t207;
t134 = t113 * t54 - t166 * t116 + t71 * t179 + t56 * t193 + t195 * t219;
t102 = sin(pkin(12));
t211 = cos(pkin(12));
t171 = t211 * t107;
t146 = -t102 * t101 + t110 * t171;
t173 = t105 * t211;
t128 = -t104 * t173 + t146 * t109;
t172 = t211 * t101;
t145 = t102 * t107 + t110 * t172;
t121 = t145 * t100 - t128 * t106;
t129 = t146 * t104 + t109 * t173;
t233 = t129 * t103 + t121 * t108;
t45 = t128 * t100 + t145 * t106;
t23 = t113 * t45 + t233 * t116;
t204 = t102 * t110;
t144 = -t107 * t204 - t172;
t131 = t102 * t199 + t144 * t109;
t143 = -t101 * t204 + t171;
t122 = t143 * t100 - t131 * t106;
t205 = t102 * t105;
t132 = t144 * t104 - t109 * t205;
t232 = t132 * t103 + t122 * t108;
t46 = t131 * t100 + t143 * t106;
t25 = t113 * t46 + t232 * t116;
t200 = t104 * t110;
t133 = t106 * t200 + t142;
t130 = t133 * t108;
t202 = t103 * t116;
t64 = t105 * t206 + (t105 * t197 + t200) * t100;
t212 = t64 * t113;
t74 = t109 * t110 - t181;
t34 = -t116 * t130 - t74 * t202 + t212;
t150 = g(1) * t25 + g(2) * t23 + g(3) * t34;
t31 = t165 * t113 + t116 * t56;
t236 = qJD(4) * t31 - t134 + t150;
t112 = sin(qJ(5));
t196 = qJD(4) * t112;
t235 = -qJD(6) * t196 + qJDD(5);
t111 = sin(qJ(6));
t114 = cos(qJ(6));
t185 = t112 * qJDD(4);
t61 = ((qJD(6) + t194) * qJD(5) + t185) * t111 - t235 * t114;
t124 = -t103 * t133 + t74 * t108;
t123 = t74 * t103 + t130;
t35 = t113 * t123 + t64 * t116;
t21 = t35 * t112 - t115 * t124;
t24 = -t233 * t113 + t45 * t116;
t26 = -t232 * t113 + t46 * t116;
t38 = t103 * t121 - t108 * t129;
t39 = t103 * t122 - t108 * t132;
t151 = g(1) * (-t112 * t26 + t115 * t39) + g(2) * (-t112 * t24 + t115 * t38) - g(3) * t21;
t167 = pkin(5) * t112 - pkin(11) * t115;
t228 = -t166 * t113 - t116 * t54;
t15 = qJDD(4) * pkin(10) - qJD(4) * t140 - t228;
t29 = qJD(4) * pkin(10) + t31;
t41 = -t103 * t55 + t108 * t71;
t20 = t112 * t41 + t115 * t29;
t40 = -t103 * t53 + t108 * t70;
t215 = t115 * t40;
t2 = -qJDD(5) * pkin(5) + qJD(5) * t20 + t112 * t15 - t215;
t231 = (pkin(11) * qJD(6) + t167 * qJD(4)) * t90 - t151 - t2;
t186 = qJD(4) * qJD(5);
t95 = t115 * qJDD(4);
t77 = t112 * t186 + qJDD(6) - t95;
t82 = t167 * qJD(5);
t84 = -pkin(5) * t115 - pkin(11) * t112 - pkin(4);
t230 = t84 * t77 + (t31 - t82) * t90;
t198 = t106 * t108;
t229 = (-t100 * t113 + t116 * t198) * t104 + t109 * t202;
t18 = qJD(5) * pkin(11) + t20;
t226 = (pkin(10) * t90 + t18) * qJD(6) + t150;
t187 = t114 * qJD(5);
t78 = t111 * t196 - t187;
t225 = t78 * t90;
t192 = qJD(5) * t111;
t80 = t114 * t196 + t192;
t224 = t80 * t90;
t98 = t112 ^ 2;
t221 = -t115 ^ 2 + t98;
t220 = qJD(4) * pkin(4);
t177 = t115 * t186;
t60 = qJD(6) * t187 + (t177 + t185) * t114 + t235 * t111;
t218 = t111 * t60;
t217 = t112 * t40;
t216 = t114 * t80;
t214 = t115 * t90;
t210 = qJD(5) * t78;
t209 = qJD(5) * t80;
t203 = t103 * t113;
t191 = qJD(5) * t112;
t190 = qJD(5) * t115;
t189 = qJD(6) * t111;
t188 = qJD(6) * t114;
t184 = t90 * t192;
t183 = t90 * t187;
t178 = t103 * t193;
t176 = t116 * t186;
t28 = t140 - t220;
t169 = -qJD(4) * t28 - t15;
t27 = t84 * qJD(4) + t140;
t4 = t111 * t27 + t114 * t18;
t164 = t111 * t18 - t114 * t27;
t22 = t112 * t124 + t35 * t115;
t10 = t111 * t34 + t114 * t22;
t9 = -t111 * t22 + t114 * t34;
t67 = t109 * t203 + (t100 * t116 + t113 * t198) * t104;
t73 = -t103 * t201 + t108 * t109;
t48 = t112 * t73 + t115 * t67;
t163 = -t111 * t229 + t114 * t48;
t162 = -t111 * t48 - t114 * t229;
t161 = t112 * t29 - t115 * t41;
t47 = t112 * t67 - t115 * t73;
t118 = qJD(4) ^ 2;
t160 = qJDD(4) * t116 - t113 * t118;
t158 = t111 * t77 - t90 * t188;
t157 = t114 * t77 + t90 * t189;
t76 = t108 * t112 + t115 * t203;
t156 = -t111 * t76 - t114 * t202;
t155 = t111 * t202 - t114 * t76;
t75 = -t108 * t115 + t112 * t203;
t149 = g(1) * t26 + g(2) * t24 + g(3) * t35;
t148 = -g(1) * t205 + g(2) * t173 - g(3) * t110;
t139 = -pkin(10) * qJDD(5) + (t28 - t140 - t220) * qJD(5);
t137 = qJD(6) * t84 * t90 - t149;
t17 = -qJD(5) * pkin(5) + t161;
t136 = -pkin(11) * t77 + (-t17 + t161) * t90;
t1 = qJDD(5) * pkin(11) - t161 * qJD(5) + t115 * t15 + t217;
t135 = -pkin(10) * t77 + qJD(5) * t17 + qJD(6) * t27 + t140 * t90 + t1;
t117 = qJD(5) ^ 2;
t127 = 0.2e1 * qJDD(4) * pkin(4) - pkin(10) * t117 + t236;
t69 = t76 * qJD(5) + t112 * t178;
t68 = -t75 * qJD(5) + t115 * t178;
t63 = t67 * qJD(4);
t62 = t229 * qJD(4);
t37 = t48 * qJD(5) + t112 * t62;
t36 = -t47 * qJD(5) + t115 * t62;
t33 = t35 * qJD(4);
t32 = (t116 * t123 - t212) * qJD(4);
t14 = t112 * t39 + t115 * t26;
t12 = t112 * t38 + t115 * t24;
t8 = -qJD(5) * t21 + t32 * t115;
t7 = qJD(5) * t22 + t32 * t112;
t6 = qJD(4) * t82 + t84 * qJDD(4) + t134;
t5 = t114 * t6;
t3 = [qJDD(1) - g(3), t110 * t86 - g(3) + (t101 ^ 2 + t107 ^ 2) * t105 ^ 2 * qJDD(1), t133 * t53 + t54 * t64 + t70 * t74 - g(3), 0, -qJD(4) * t33 - qJDD(4) * t34, -qJD(4) * t32 - qJDD(4) * t35, 0, 0, 0, 0, 0, -t34 * t95 - qJD(5) * t7 - qJDD(5) * t21 + (-t115 * t33 + t34 * t191) * qJD(4), t34 * t185 - qJD(5) * t8 - qJDD(5) * t22 + (t112 * t33 + t34 * t190) * qJD(4), 0, 0, 0, 0, 0 -(-qJD(6) * t10 - t111 * t8 + t114 * t33) * t90 + t9 * t77 + t7 * t78 + t21 * t61 (qJD(6) * t9 + t111 * t33 + t114 * t8) * t90 - t10 * t77 + t7 * t80 + t21 * t60; 0, t148 + t86, t109 * t70 + (t100 * t54 + t106 * t53) * t104 + t148, 0, -qJD(4) * t63 + qJDD(4) * t229, -qJD(4) * t62 - qJDD(4) * t67, 0, 0, 0, 0, 0, t229 * t95 - qJD(5) * t37 - qJDD(5) * t47 + (-t115 * t63 - t191 * t229) * qJD(4), -t229 * t185 - qJD(5) * t36 - qJDD(5) * t48 + (t112 * t63 - t190 * t229) * qJD(4), 0, 0, 0, 0, 0 -(-qJD(6) * t163 - t111 * t36 + t114 * t63) * t90 + t162 * t77 + t37 * t78 + t47 * t61 (qJD(6) * t162 + t111 * t63 + t114 * t36) * t90 - t163 * t77 + t37 * t80 + t47 * t60; 0, 0, g(1) * t132 + g(2) * t129 - g(3) * t74 + t70, 0, t160 * t103 (-qJDD(4) * t113 - t116 * t118) * t103, 0, 0, 0, 0, 0, -qJD(5) * t69 - qJDD(5) * t75 + (-t112 * t176 + t115 * t160) * t103, -qJD(5) * t68 - qJDD(5) * t76 + (-t112 * t160 - t115 * t176) * t103, 0, 0, 0, 0, 0 -(qJD(6) * t155 - t111 * t68 + t114 * t179) * t90 + t156 * t77 + t69 * t78 + t75 * t61 (qJD(6) * t156 + t111 * t179 + t114 * t68) * t90 + t155 * t77 + t69 * t80 + t75 * t60; 0, 0, 0, qJDD(4), t236, t149 + t228, qJDD(4) * t98 + 0.2e1 * t112 * t177, 0.2e1 * t112 * t95 - 0.2e1 * t221 * t186, qJDD(5) * t112 + t115 * t117, qJDD(5) * t115 - t112 * t117, 0, t112 * t139 + t115 * t127, -t112 * t127 + t115 * t139, t115 * t80 * t187 + (t114 * t60 - t189 * t80) * t112 (-t111 * t80 - t114 * t78) * t190 + (-t218 - t114 * t61 + (t111 * t78 - t216) * qJD(6)) * t112 (-t60 - t183) * t115 + (t157 + t209) * t112 (t61 + t184) * t115 + (-t158 - t210) * t112, -t115 * t77 - t191 * t90, t230 * t114 + t137 * t111 + (pkin(10) * t210 + t135 * t111 + t226 * t114 - t5) * t115 + (t17 * t188 - t164 * qJD(5) + t2 * t111 + t140 * t78 + (t61 - t184) * pkin(10)) * t112, -t230 * t111 + t137 * t114 + (pkin(10) * t209 + t135 * t114 + (-t226 + t6) * t111) * t115 + (-t17 * t189 - t4 * qJD(5) + t2 * t114 + t140 * t80 + (t60 - t183) * pkin(10)) * t112; 0, 0, 0, 0, 0, 0, -t112 * t118 * t115, t221 * t118, t185, t95, qJDD(5), t112 * t169 - t151 + t215, g(1) * t14 + g(2) * t12 + g(3) * t22 + t169 * t115 - t217, -t90 * t216 + t218 (t60 + t225) * t114 + (-t61 + t224) * t111 (-t112 * t80 + t114 * t214) * qJD(4) + t158 (-t111 * t214 + t112 * t78) * qJD(4) + t157, t90 * t196, -pkin(5) * t61 + t136 * t111 + t231 * t114 + t164 * t196 - t20 * t78, -pkin(5) * t60 - t231 * t111 + t136 * t114 + t4 * t196 - t20 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 * t80, -t78 ^ 2 + t80 ^ 2, t60 - t225, -t224 - t61, t77, -t111 * t1 + t5 - t17 * t80 - g(1) * (-t111 * t14 + t114 * t25) - g(2) * (-t111 * t12 + t114 * t23) - g(3) * t9 - t237 * t4, -t114 * t1 - t111 * t6 + t17 * t78 - g(1) * (-t111 * t25 - t114 * t14) - g(2) * (-t111 * t23 - t114 * t12) + g(3) * t10 + t237 * t164;];
tau_reg  = t3;
