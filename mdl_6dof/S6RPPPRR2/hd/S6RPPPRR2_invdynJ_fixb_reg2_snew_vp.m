% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPPRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:36:16
% EndTime: 2019-05-05 13:36:24
% DurationCPUTime: 2.67s
% Computational Cost: add. (7960->265), mult. (16987->376), div. (0->0), fcn. (11447->10), ass. (0->174)
t148 = sin(qJ(5));
t141 = sin(pkin(10));
t151 = cos(qJ(5));
t143 = cos(pkin(10));
t191 = t143 * t148;
t165 = t141 * t151 + t191;
t110 = t165 * qJD(1);
t112 = (-t141 * t148 + t143 * t151) * qJD(1);
t193 = t112 * t110;
t216 = qJDD(5) - t193;
t219 = t148 * t216;
t218 = t151 * t216;
t140 = qJDD(1) * pkin(2);
t154 = qJD(1) ^ 2;
t149 = sin(qJ(1));
t152 = cos(qJ(1));
t173 = t149 * g(1) - g(2) * t152;
t118 = qJDD(1) * pkin(1) + t173;
t166 = g(1) * t152 + g(2) * t149;
t119 = -pkin(1) * t154 - t166;
t142 = sin(pkin(9));
t144 = cos(pkin(9));
t170 = -t144 * t118 + t119 * t142;
t84 = -t154 * qJ(3) + qJDD(3) - t140 + t170;
t159 = -0.2e1 * qJD(4) * qJD(1) + t84;
t217 = -qJDD(1) * qJ(4) + t159;
t147 = sin(qJ(6));
t150 = cos(qJ(6));
t100 = qJD(5) * t147 + t112 * t150;
t98 = -t150 * qJD(5) + t112 * t147;
t77 = t100 * t98;
t186 = t112 * qJD(5);
t209 = t165 * qJDD(1);
t92 = -t209 - t186;
t85 = qJDD(6) - t92;
t212 = -t77 + t85;
t215 = t147 * t212;
t214 = t150 * t212;
t190 = -g(3) + qJDD(2);
t171 = t141 * t190;
t205 = t141 * pkin(4);
t157 = -t171 + (-t154 * t205 + (-pkin(7) - qJ(4)) * qJDD(1) + t159) * t143;
t181 = t141 * qJDD(1);
t135 = t141 ^ 2;
t192 = t135 * t154;
t73 = t217 * t141 + t143 * t190;
t64 = -pkin(4) * t192 - pkin(7) * t181 + t73;
t39 = t148 * t157 + t151 * t64;
t168 = pkin(1) * t144 + pkin(2) + qJ(4);
t213 = -pkin(7) - t168;
t136 = t143 ^ 2;
t188 = t135 + t136;
t210 = t188 * t154;
t208 = -t154 * qJ(4) + qJDD(4);
t106 = qJD(6) + t110;
t180 = t143 * qJDD(1);
t109 = -t148 * t181 + t151 * t180;
t187 = qJD(5) * t110;
t94 = t109 - t187;
t172 = -t150 * qJDD(5) + t147 * t94;
t48 = (qJD(6) - t106) * t100 + t172;
t96 = t98 ^ 2;
t97 = t100 ^ 2;
t105 = t106 ^ 2;
t107 = t110 ^ 2;
t108 = t112 ^ 2;
t206 = pkin(1) * t142;
t153 = qJD(5) ^ 2;
t38 = t148 * t64 - t151 * t157;
t86 = pkin(5) * t110 - pkin(8) * t112;
t29 = -qJDD(5) * pkin(5) - t153 * pkin(8) + t112 * t86 + t38;
t204 = t147 * t29;
t56 = t77 + t85;
t203 = t147 * t56;
t132 = 0.2e1 * qJD(3) * qJD(1);
t189 = t142 * t118 + t144 * t119;
t169 = -t154 * pkin(2) + t132 + t189;
t174 = qJ(3) + t205;
t68 = t174 * qJDD(1) + t169 + t208 + (-t136 * t154 - t192) * pkin(7);
t202 = t148 * t68;
t89 = qJDD(5) + t193;
t201 = t148 * t89;
t200 = t150 * t29;
t199 = t150 * t56;
t197 = t151 * t68;
t196 = t151 * t89;
t195 = t106 * t147;
t194 = t106 * t150;
t184 = qJD(6) + t106;
t183 = qJDD(1) * qJ(3);
t179 = t148 * t77;
t178 = t151 * t77;
t176 = -pkin(5) * t151 - pkin(4);
t175 = qJ(3) + t206;
t30 = -t153 * pkin(5) + qJDD(5) * pkin(8) - t110 * t86 + t39;
t32 = (-t94 + t187) * pkin(8) + (-t92 + t186) * pkin(5) + t68;
t14 = t147 * t30 - t150 * t32;
t15 = t147 * t32 + t150 * t30;
t6 = t14 * t147 + t150 * t15;
t19 = t148 * t38 + t151 * t39;
t167 = -pkin(1) * (qJDD(1) * t142 + t144 * t154) - t189;
t3 = t148 * t6 - t151 * t29;
t4 = t148 * t29 + t151 * t6;
t1 = t141 * t4 + t143 * t3;
t5 = -t14 * t150 + t147 * t15;
t18 = t148 * t39 - t151 * t38;
t7 = t141 * t19 + t143 * t18;
t72 = -t217 * t143 + t171;
t45 = t141 * t73 - t143 * t72;
t164 = -qJDD(5) * t147 - t150 * t94;
t163 = t174 + t206;
t162 = -pkin(1) * (-qJDD(1) * t144 + t142 * t154) - t170;
t83 = t169 + t183;
t79 = t83 + t208;
t161 = t175 * qJDD(1) + t79;
t66 = -qJD(6) * t98 - t164;
t120 = t188 * qJDD(1);
t117 = t141 * t210;
t116 = t143 * t210;
t103 = -t108 - t153;
t102 = -t108 + t153;
t101 = t107 - t153;
t93 = t109 - 0.2e1 * t187;
t91 = t209 + 0.2e1 * t186;
t87 = -t153 - t107;
t82 = t106 * t98;
t81 = -t97 + t105;
t80 = t96 - t105;
t76 = -t107 - t108;
t75 = t97 - t96;
t71 = -t97 - t105;
t70 = -t103 * t148 - t196;
t69 = t103 * t151 - t201;
t67 = -t105 - t96;
t65 = -qJD(6) * t100 - t172;
t63 = t96 + t97;
t62 = t109 * t148 - t151 * t209;
t61 = -t109 * t151 - t148 * t209;
t60 = t151 * t87 - t219;
t59 = t148 * t87 + t218;
t54 = (t100 * t147 - t150 * t98) * t106;
t53 = t184 * t98 + t164;
t52 = t66 + t82;
t51 = t66 - t82;
t49 = -t184 * t100 - t172;
t47 = -t100 * t195 + t150 * t66;
t46 = -t147 * t65 + t98 * t194;
t44 = t141 * t70 + t143 * t69;
t43 = t150 * t80 - t203;
t42 = -t147 * t81 + t214;
t41 = -t147 * t71 - t199;
t40 = t150 * t71 - t203;
t36 = t141 * t62 + t143 * t61;
t35 = t150 * t67 - t215;
t34 = t147 * t67 + t214;
t33 = t141 * t60 + t143 * t59;
t28 = t147 * t52 - t150 * t48;
t27 = -t147 * t51 + t150 * t49;
t26 = -t147 * t48 - t150 * t52;
t25 = -t148 * t53 + t151 * t41;
t24 = t148 * t41 + t151 * t53;
t23 = -t148 * t49 + t151 * t35;
t22 = t148 * t35 + t151 * t49;
t21 = -t148 * t63 + t151 * t28;
t20 = t148 * t28 + t151 * t63;
t17 = -pkin(8) * t40 + t200;
t16 = -pkin(8) * t34 + t204;
t12 = -pkin(5) * t40 + t15;
t11 = t141 * t25 + t143 * t24;
t10 = -pkin(5) * t34 + t14;
t9 = t141 * t23 + t143 * t22;
t8 = t141 * t21 + t143 * t20;
t2 = -pkin(8) * t26 - t5;
t13 = [0, 0, 0, 0, 0, qJDD(1), t173, t166, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t162, t167, 0, pkin(1) * (t142 * t189 - t144 * t170), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t140 - t162, t132 - t167 + 0.2e1 * t183, pkin(1) * (t142 * t83 - t144 * t84) - pkin(2) * t84 + qJ(3) * t83, t136 * qJDD(1), -0.2e1 * t141 * t180, 0, t135 * qJDD(1), 0, 0, t168 * t117 + t161 * t141, t168 * t116 + t161 * t143, t168 * t120 - t175 * t210 - t45, -t168 * t45 + t175 * t79, -t141 * (t148 * t94 + t151 * t186) + t143 * (-t148 * t186 + t151 * t94), -t141 * (-t148 * t91 + t151 * t93) + t143 * (-t148 * t93 - t151 * t91), -t141 * (t102 * t151 + t219) + t143 * (-t102 * t148 + t218), -t141 * (t148 * t187 + t151 * t92) + t143 * (-t148 * t92 + t151 * t187), -t141 * (t101 * t148 + t196) + t143 * (t101 * t151 - t201), (-t141 * (-t110 * t148 - t112 * t151) + t143 * (-t110 * t151 + t112 * t148)) * qJD(5), -t141 * (pkin(7) * t60 - t197) + t143 * (-pkin(7) * t59 + t202) + t163 * t91 - t168 * t33, -t141 * (pkin(7) * t70 + t202) + t143 * (-pkin(7) * t69 + t197) + t163 * t93 - t168 * t44, -t141 * (pkin(7) * t62 + t19) + t143 * (-pkin(7) * t61 - t18) + t163 * t76 - t168 * t36, t163 * t68 + t213 * t7, -t141 * (t148 * t47 - t178) + t143 * (t151 * t47 + t179), -t141 * (t148 * t27 - t151 * t75) + t143 * (t148 * t75 + t151 * t27), -t141 * (t148 * t42 - t151 * t52) + t143 * (t148 * t52 + t151 * t42), -t141 * (t148 * t46 + t178) + t143 * (t151 * t46 - t179), -t141 * (t148 * t43 + t151 * t48) + t143 * (-t148 * t48 + t151 * t43), -t141 * (t148 * t54 - t151 * t85) + t143 * (t148 * t85 + t151 * t54), -t141 * (pkin(7) * t23 + t10 * t151 + t148 * t16) + t143 * (-pkin(7) * t22 - t10 * t148 + t151 * t16) + t163 * t34 - t168 * t9, -t141 * (pkin(7) * t25 + t12 * t151 + t148 * t17) + t143 * (-pkin(7) * t24 - t12 * t148 + t151 * t17) + t163 * t40 - t168 * t11, -t141 * (pkin(7) * t21 + t148 * t2) + t143 * (-pkin(7) * t20 + t151 * t2) - t168 * t8 + (pkin(5) * t191 - t141 * t176 + t175) * t26, (-t141 * (-pkin(8) * t148 + t176) + t143 * (pkin(5) * t148 - pkin(8) * t151) + t175) * t5 + t213 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141 * t72 + t143 * t73, 0, 0, 0, 0, 0, 0, -t141 * t59 + t143 * t60, -t141 * t69 + t143 * t70, -t141 * t61 + t143 * t62, -t141 * t18 + t143 * t19, 0, 0, 0, 0, 0, 0, -t141 * t22 + t143 * t23, -t141 * t24 + t143 * t25, -t141 * t20 + t143 * t21, -t141 * t3 + t143 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t154, t84, 0, 0, 0, 0, 0, 0, -t117, -t116, -t120, t45, 0, 0, 0, 0, 0, 0, t33, t44, t36, t7, 0, 0, 0, 0, 0, 0, t9, t11, t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, t180, -t210, t79, 0, 0, 0, 0, 0, 0, t91, t93, t76, t68, 0, 0, 0, 0, 0, 0, t34, t40, t26, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t108 - t107, t109, -t193, -t209, qJDD(5), -t38, -t39, 0, 0, t100 * t194 + t147 * t66, t147 * t49 + t150 * t51, t150 * t81 + t215, t150 * t65 + t195 * t98, t147 * t80 + t199, (-t100 * t150 - t147 * t98) * t106, pkin(5) * t49 + pkin(8) * t35 - t200, pkin(5) * t53 + pkin(8) * t41 + t204, pkin(5) * t63 + pkin(8) * t28 + t6, -pkin(5) * t29 + pkin(8) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t75, t52, -t77, -t48, t85, -t14, -t15, 0, 0;];
tauJ_reg  = t13;
