% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:40
% EndTime: 2019-03-09 01:51:43
% DurationCPUTime: 1.68s
% Computational Cost: add. (1127->294), mult. (1925->340), div. (0->0), fcn. (1035->6), ass. (0->171)
t92 = pkin(1) + qJ(3);
t204 = qJD(1) * t92;
t53 = -qJD(2) + t204;
t167 = t53 * qJD(1);
t98 = cos(qJ(1));
t82 = g(1) * t98;
t95 = sin(qJ(1));
t179 = g(2) * t95 + t82;
t212 = t179 + t167;
t68 = qJ(2) * qJD(1) + qJD(3);
t52 = -pkin(7) * qJD(1) + t68;
t97 = cos(qJ(4));
t22 = (pkin(5) * qJD(1) - t52) * t97;
t162 = qJD(5) + t22;
t151 = qJDD(4) * qJ(5);
t85 = qJD(1) * qJD(2);
t86 = qJ(2) * qJDD(1);
t141 = qJDD(3) + t85 + t86;
t40 = -pkin(7) * qJDD(1) + t141;
t94 = sin(qJ(4));
t28 = t94 * t40;
t128 = -t28 - t151;
t185 = t97 * t52;
t10 = (-qJD(5) - t185) * qJD(4) + t128;
t129 = -t97 * t40 + qJDD(5);
t172 = qJD(4) * t94;
t121 = t52 * t172 + t129;
t169 = qJDD(4) * pkin(4);
t11 = t121 - t169;
t131 = qJD(4) * pkin(4) - qJD(5);
t24 = -t131 - t185;
t160 = qJD(4) * qJ(5);
t45 = t94 * t52;
t27 = -t45 - t160;
t105 = (t24 * t94 - t27 * t97) * qJD(4) - t10 * t94 - t11 * t97;
t211 = t105 - t179;
t152 = qJD(1) * qJD(4);
t139 = t97 * t152;
t155 = t94 * qJDD(1);
t210 = t139 + t155;
t93 = sin(qJ(6));
t166 = t93 * qJD(4);
t175 = qJD(1) * t94;
t96 = cos(qJ(6));
t209 = t96 * t175 - t166;
t180 = g(1) * t95 - g(2) * t98;
t176 = t97 * qJ(5);
t124 = t176 - t92;
t77 = t94 * pkin(4);
t208 = t124 - t77;
t140 = t94 * t152;
t70 = t97 * qJDD(1);
t207 = -t140 + t70;
t164 = t97 * qJD(1);
t59 = qJD(6) + t164;
t191 = t59 * t93;
t35 = -qJDD(6) - t207;
t25 = t96 * t35;
t109 = -qJD(6) * t191 - t25;
t174 = qJD(4) * t209;
t206 = -t174 + t25;
t205 = -t77 + t176;
t21 = -pkin(5) * t175 + t45;
t16 = t21 + t160;
t99 = -pkin(4) - pkin(8);
t203 = -t99 * t35 + (t16 - t21) * t59;
t165 = t96 * qJD(4);
t39 = t93 * t175 + t165;
t9 = t39 * qJD(6) + t93 * qJDD(4) - t210 * t96;
t91 = -pkin(7) + qJ(2);
t159 = qJDD(4) * t91;
t202 = (qJD(2) + t53 + t204) * qJD(4) + t159;
t163 = pkin(4) * t175 - qJD(2);
t17 = -t124 * qJD(1) + t163;
t201 = (-qJD(1) * t208 + qJD(2) + t17) * qJD(4) + t159;
t75 = 0.2e1 * t85;
t200 = g(3) * t94;
t199 = g(3) * t97;
t8 = qJD(6) * t209 + t96 * qJDD(4) + t210 * t93;
t198 = t8 * t96;
t197 = pkin(5) - t91;
t196 = t27 * t94;
t195 = t209 * t59;
t194 = t209 * t94;
t193 = t39 * t59;
t192 = t39 * t94;
t190 = t93 * t35;
t189 = t93 * t94;
t188 = t94 * t96;
t187 = t95 * t97;
t186 = t96 * t97;
t184 = t98 * t93;
t183 = t98 * t96;
t182 = -g(2) * t187 - t97 * t82;
t42 = pkin(4) * t164 + qJ(5) * t175;
t181 = t98 * pkin(1) + t95 * qJ(2);
t88 = t94 ^ 2;
t89 = t97 ^ 2;
t178 = -t88 - t89;
t101 = qJD(1) ^ 2;
t177 = t88 * t101;
t173 = qJD(4) * t39;
t171 = qJD(4) * t97;
t170 = qJD(6) * t96;
t168 = t17 * qJD(1);
t100 = qJD(4) ^ 2;
t161 = t100 + t101;
t158 = qJDD(4) * t94;
t157 = qJDD(4) * t97;
t156 = t92 * qJDD(1);
t90 = qJDD(1) * pkin(1);
t154 = t90 - qJDD(2);
t153 = qJ(5) * qJDD(1);
t84 = qJD(3) * qJD(1);
t150 = t97 * t191;
t149 = t59 * t186;
t148 = t98 * qJ(3) + t181;
t147 = t97 * t101 * t94;
t144 = pkin(4) * t171 + t94 * t160 + qJD(3);
t83 = qJDD(1) * qJ(3);
t41 = t83 + t84 + t154;
t107 = pkin(4) * t210 + qJ(5) * t140 + t41;
t130 = qJD(4) * pkin(8) - qJD(5);
t3 = pkin(8) * t155 + (t130 * qJD(1) - t153) * t97 + t107;
t5 = pkin(5) * t207 + t99 * qJDD(4) + t121;
t143 = -t93 * t3 + t96 * t5;
t142 = qJDD(2) - t180;
t138 = qJD(4) * t197;
t136 = qJD(1) * t42 - g(3);
t108 = t94 * pkin(8) - t124;
t12 = t108 * qJD(1) + t163;
t135 = qJD(6) * t12 - t5;
t13 = t99 * qJD(4) + t162;
t134 = qJD(6) * t13 + t3;
t133 = 0.2e1 * t86 + t75 - t179;
t132 = t178 * qJDD(1);
t127 = qJD(6) * t97 + qJD(1);
t126 = -t90 + t142;
t125 = t182 + t200;
t120 = pkin(4) * t97 + qJ(5) * t94;
t2 = t96 * t12 + t93 * t13;
t114 = -t83 + t126;
t113 = t173 - t190;
t112 = t59 ^ 2;
t110 = t59 * t170 - t190;
t106 = -t100 * t91 + t180;
t20 = -t97 * qJD(5) + t144;
t7 = (-qJD(1) * qJD(5) - t153) * t97 + t107;
t104 = -qJD(1) * t20 + qJDD(1) * t208 - t106 - t7;
t103 = t106 + t41 + t84 + t156;
t6 = -pkin(5) * t155 + (qJD(5) - t22) * qJD(4) - t128;
t102 = -t199 + t6 - t179 * t94 + (pkin(8) * t164 - qJD(6) * t99 + t42) * t59;
t74 = t98 * qJ(2);
t71 = t89 * t101;
t48 = t197 * t97;
t47 = t197 * t94;
t44 = 0.2e1 * t139 + t155;
t43 = -t70 + 0.2e1 * t140;
t34 = -t161 * t94 + t157;
t33 = t161 * t97 + t158;
t32 = t93 * t187 - t183;
t31 = t95 * t186 + t184;
t30 = t97 * t184 + t95 * t96;
t29 = -t97 * t183 + t95 * t93;
t26 = t108 + t77;
t19 = t94 * qJD(2) - t97 * t138;
t18 = -t97 * qJD(2) - t94 * t138;
t15 = t130 * t97 + t144;
t14 = t17 * t164;
t1 = -t93 * t12 + t96 * t13;
t4 = [qJDD(1), t180, t179, -0.2e1 * t90 + t142, t133, t154 * pkin(1) - g(1) * (-t95 * pkin(1) + t74) - g(2) * t181 + (t75 + t86) * qJ(2), qJDD(3) + t133, -t114 + 0.2e1 * t84 + t156, t41 * t92 + t53 * qJD(3) + t141 * qJ(2) + t68 * qJD(2) - g(1) * (-t92 * t95 + t74) - g(2) * t148, t89 * qJDD(1) - 0.2e1 * t94 * t139, -0.2e1 * t94 * t70 + 0.2e1 * (t88 - t89) * t152, -t100 * t94 + t157, -t100 * t97 - t158, 0, t103 * t94 + t202 * t97, t103 * t97 - t202 * t94, t91 * t132 + t178 * t85 - t211, t104 * t94 - t201 * t97, t104 * t97 + t201 * t94, -t7 * t208 + t17 * t20 - g(1) * (-t98 * pkin(7) + t74) - g(2) * (-t205 * t98 + t148) + (-t24 * t97 - t196) * qJD(2) + (g(2) * pkin(7) - g(1) * t208) * t95 + t105 * t91, t8 * t189 + (t97 * t166 + t94 * t170) * t39 (t209 * t93 + t39 * t96) * t171 + (t198 - t9 * t93 + (t209 * t96 - t39 * t93) * qJD(6)) * t94 (t59 * t166 + t8) * t97 + (t110 - t173) * t94 (t59 * t165 - t9) * t97 + (t109 - t174) * t94, -t59 * t172 - t35 * t97 (-t93 * t15 + t96 * t18) * t59 - (-t93 * t26 + t96 * t48) * t35 + t143 * t97 - t19 * t209 - t47 * t9 - t6 * t188 - g(1) * t32 + g(2) * t30 + (-t1 * t94 - t16 * t186) * qJD(4) + ((-t96 * t26 - t93 * t48) * t59 - t2 * t97 + t16 * t189) * qJD(6), t2 * t172 - g(1) * t31 - g(2) * t29 + t19 * t39 - t47 * t8 + (-(qJD(6) * t48 + t15) * t59 + t26 * t35 - t134 * t97 + t16 * qJD(6) * t94) * t96 + (-(-qJD(6) * t26 + t18) * t59 + t48 * t35 + t6 * t94 + (t16 * qJD(4) + t135) * t97) * t93; 0, 0, 0, qJDD(1), -t101, -t101 * qJ(2) + t126, -t101, -qJDD(1), -t68 * qJD(1) + t114 - t84, 0, 0, 0, 0, 0, -t44, t43, t71 + t177, t44, -t43, t97 * t153 + (t196 + (qJD(5) + t24) * t97) * qJD(1) - t107 - t180, 0, 0, 0, 0, 0 (t149 + t194) * qJD(1) + t110 (-t150 - t192) * qJD(1) + t109; 0, 0, 0, 0, 0, 0, qJDD(1), -t101, t141 - t212, 0, 0, 0, 0, 0, t34, -t33, t132, -t34, t33, -t168 + t211, 0, 0, 0, 0, 0, t94 * t9 + t206 * t97 + (t127 * t93 + t165 * t94) * t59, t94 * t8 + t113 * t97 + (t127 * t96 - t166 * t94) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t71 - t177, t70, -t155, qJDD(4) (t40 - t167) * t97 + t125, t212 * t94 + t199 - t28, -t120 * qJDD(1) + ((-t27 - t160) * t97 + (t131 + t24) * t94) * qJD(1), t136 * t94 + t129 + t14 - 0.2e1 * t169 - t182, 0.2e1 * t151 + 0.2e1 * qJD(4) * qJD(5) + t28 + t136 * t97 + (-t179 - t168) * t94, -t10 * qJ(5) - t11 * pkin(4) - t17 * t42 - t24 * t45 - g(3) * t205 + (-qJD(5) + t185) * t27 - t179 * t120, -t193 * t93 + t198 (-t9 - t193) * t96 + (-t8 - t195) * t93 (-t150 + t192) * qJD(1) + t109 (-t149 + t194) * qJD(1) - t110, t59 * t175, qJ(5) * t9 + t1 * t175 + t102 * t93 - t162 * t209 + t203 * t96, qJ(5) * t8 + t102 * t96 + t162 * t39 - t2 * t175 - t203 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, qJDD(4) - t147, -t71 - t100, t27 * qJD(4) + t11 - t125 + t14, 0, 0, 0, 0, 0, -t112 * t93 - t206, -t112 * t96 - t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39 * t209, -t209 ^ 2 + t39 ^ 2, t8 - t195, t193 - t9, -t35, -g(1) * t29 + g(2) * t31 - g(3) * t188 - t16 * t39 + t143 + (-qJD(6) + t59) * t2, -g(1) * t30 - g(2) * t32 + t1 * t59 - t16 * t209 - t134 * t96 + (t135 + t200) * t93;];
tau_reg  = t4;
