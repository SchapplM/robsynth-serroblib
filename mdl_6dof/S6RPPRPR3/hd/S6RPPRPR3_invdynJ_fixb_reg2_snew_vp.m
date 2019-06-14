% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:09:51
% EndTime: 2019-05-05 14:09:59
% DurationCPUTime: 3.00s
% Computational Cost: add. (9197->292), mult. (19153->411), div. (0->0), fcn. (12338->10), ass. (0->187)
t152 = sin(pkin(10));
t154 = cos(pkin(10));
t160 = sin(qJ(4));
t163 = cos(qJ(4));
t203 = t163 * t152;
t115 = (t160 * t154 + t203) * qJD(1);
t199 = qJD(1) * t163;
t117 = -t152 * t160 * qJD(1) + t154 * t199;
t208 = t117 * t115;
t225 = qJDD(4) - t208;
t227 = t152 * t225;
t226 = t154 * t225;
t198 = qJD(4) * t115;
t191 = qJD(1) * qJD(4);
t183 = t163 * t191;
t189 = t160 * qJDD(1);
t126 = -t183 - t189;
t141 = t163 * qJDD(1);
t184 = t160 * t191;
t127 = t141 - t184;
t97 = t152 * t126 + t154 * t127;
t81 = t97 - t198;
t159 = sin(qJ(6));
t162 = cos(qJ(6));
t100 = -t162 * qJD(4) + t159 * t117;
t102 = t159 * qJD(4) + t162 * t117;
t78 = t102 * t100;
t177 = -t154 * t126 + t152 * t127;
t94 = qJDD(6) + t177;
t222 = -t78 + t94;
t224 = t159 * t222;
t223 = t162 * t222;
t197 = qJD(4) * t117;
t79 = t177 + t197;
t155 = cos(pkin(9));
t176 = pkin(1) * t155 + pkin(2) + pkin(7);
t221 = -qJ(5) - t176;
t149 = -g(3) + qJDD(2);
t151 = qJDD(1) * pkin(2);
t166 = qJD(1) ^ 2;
t161 = sin(qJ(1));
t164 = cos(qJ(1));
t181 = t161 * g(1) - t164 * g(2);
t123 = qJDD(1) * pkin(1) + t181;
t174 = t164 * g(1) + t161 * g(2);
t124 = -t166 * pkin(1) - t174;
t153 = sin(pkin(9));
t178 = t155 * t123 - t153 * t124;
t88 = -t166 * qJ(3) + qJDD(3) - t151 - t178;
t169 = -qJDD(1) * pkin(7) + t88;
t168 = t163 * t169;
t202 = t163 * t166;
t167 = t168 - t127 * qJ(5) + qJDD(4) * pkin(4) + (-pkin(4) * t202 - qJ(5) * t191 - t149) * t160;
t133 = qJD(4) * pkin(4) - qJ(5) * t199;
t147 = t160 ^ 2;
t207 = t147 * t166;
t75 = t163 * t149 + t160 * t169;
t60 = -pkin(4) * t207 + t126 * qJ(5) - qJD(4) * t133 + t75;
t34 = -0.2e1 * qJD(5) * t115 + t152 * t167 + t154 * t60;
t112 = qJD(6) + t115;
t179 = -t162 * qJDD(4) + t159 * t97;
t48 = (qJD(6) - t112) * t102 + t179;
t98 = t100 ^ 2;
t99 = t102 ^ 2;
t111 = t112 ^ 2;
t113 = t115 ^ 2;
t114 = t117 ^ 2;
t220 = 0.2e1 * qJD(5);
t143 = 0.2e1 * qJD(3) * qJD(1);
t190 = qJDD(1) * qJ(3);
t201 = t153 * t123 + t155 * t124;
t87 = -t166 * pkin(2) + t143 + t190 + t201;
t83 = -t166 * pkin(7) + t87;
t64 = -t126 * pkin(4) - qJ(5) * t207 + t133 * t199 + qJDD(5) + t83;
t218 = t152 * t64;
t92 = qJDD(4) + t208;
t217 = t152 * t92;
t216 = t154 * t64;
t215 = t154 * t92;
t165 = qJD(4) ^ 2;
t180 = t152 * t60 - t154 * t167;
t89 = t115 * pkin(5) - t117 * pkin(8);
t26 = -qJDD(4) * pkin(5) - t165 * pkin(8) + (t220 + t89) * t117 + t180;
t214 = t159 * t26;
t62 = t78 + t94;
t213 = t159 * t62;
t212 = t162 * t26;
t211 = t162 * t62;
t210 = t112 * t159;
t209 = t112 * t162;
t148 = t163 ^ 2;
t206 = t148 * t166;
t186 = t160 * t202;
t134 = qJDD(4) + t186;
t205 = t160 * t134;
t135 = qJDD(4) - t186;
t204 = t163 * t135;
t200 = t147 + t148;
t196 = qJD(4) * t152;
t195 = qJD(4) * t154;
t192 = qJD(6) + t112;
t188 = t152 * t78;
t187 = t154 * t78;
t185 = -pkin(5) * t154 - pkin(4);
t182 = pkin(1) * t153 + qJ(3);
t27 = -t165 * pkin(5) + qJDD(4) * pkin(8) - t115 * t89 + t34;
t36 = t79 * pkin(5) - t81 * pkin(8) + t64;
t14 = t159 * t27 - t162 * t36;
t15 = t159 * t36 + t162 * t27;
t6 = t159 * t14 + t162 * t15;
t33 = t117 * t220 + t180;
t19 = t152 * t33 + t154 * t34;
t175 = -pkin(1) * (t153 * qJDD(1) + t155 * t166) - t201;
t2 = t152 * t6 - t154 * t26;
t3 = t152 * t26 + t154 * t6;
t1 = t160 * t3 + t163 * t2;
t5 = -t162 * t14 + t159 * t15;
t18 = t152 * t34 - t154 * t33;
t7 = t160 * t19 + t163 * t18;
t74 = t160 * t149 - t168;
t47 = t160 * t75 - t163 * t74;
t172 = -t159 * qJDD(4) - t162 * t97;
t171 = t160 * pkin(4) + t182;
t170 = -pkin(1) * (-t155 * qJDD(1) + t153 * t166) + t178;
t80 = -t177 + t197;
t70 = -t100 * qJD(6) - t172;
t138 = -t165 - t206;
t137 = -t165 - t207;
t131 = t200 * qJDD(1);
t128 = t141 - 0.2e1 * t184;
t125 = 0.2e1 * t183 + t189;
t107 = -t114 - t165;
t106 = -t114 + t165;
t105 = t113 - t165;
t104 = t163 * t138 - t205;
t103 = t160 * t137 + t204;
t90 = -t165 - t113;
t86 = t112 * t100;
t85 = -t99 + t111;
t84 = t98 - t111;
t82 = t198 + t97;
t77 = -t113 - t114;
t76 = t99 - t98;
t73 = -t99 - t111;
t72 = -t152 * t107 - t215;
t71 = t154 * t107 - t217;
t69 = -t102 * qJD(6) - t179;
t68 = -t111 - t98;
t67 = t98 + t99;
t66 = t154 * t90 - t227;
t65 = t152 * t90 + t226;
t56 = (-t100 * t162 + t102 * t159) * t112;
t55 = t152 * t82 + t154 * t80;
t54 = t152 * t80 - t154 * t82;
t53 = t192 * t100 + t172;
t52 = t70 + t86;
t51 = t70 - t86;
t49 = -t192 * t102 - t179;
t46 = -t102 * t210 + t162 * t70;
t45 = t100 * t209 - t159 * t69;
t44 = t160 * t72 + t163 * t71;
t43 = t162 * t84 - t213;
t42 = -t159 * t85 + t223;
t41 = -t159 * t73 - t211;
t40 = t162 * t73 - t213;
t39 = t162 * t68 - t224;
t38 = t159 * t68 + t223;
t37 = t160 * t66 + t163 * t65;
t31 = t160 * t55 + t163 * t54;
t30 = t159 * t52 - t162 * t48;
t29 = -t159 * t51 + t162 * t49;
t28 = -t159 * t48 - t162 * t52;
t25 = -t152 * t53 + t154 * t41;
t24 = t152 * t41 + t154 * t53;
t23 = -t152 * t49 + t154 * t39;
t22 = t152 * t39 + t154 * t49;
t21 = -t152 * t67 + t154 * t30;
t20 = t152 * t30 + t154 * t67;
t17 = -pkin(8) * t40 + t212;
t16 = -pkin(8) * t38 + t214;
t12 = t160 * t25 + t163 * t24;
t11 = t160 * t23 + t163 * t22;
t10 = -pkin(5) * t40 + t15;
t9 = -pkin(5) * t38 + t14;
t8 = t160 * t21 + t163 * t20;
t4 = -pkin(8) * t28 - t5;
t13 = [0, 0, 0, 0, 0, qJDD(1), t181, t174, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t170, t175, 0, pkin(1) * (t153 * t201 + t155 * t178), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t151 - t170, t143 - t175 + 0.2e1 * t190, pkin(1) * (t153 * t87 - t155 * t88) - pkin(2) * t88 + qJ(3) * t87, (t127 - t184) * t163, -t163 * t125 - t160 * t128, t204 - t160 * (t165 - t206), (-t126 + t183) * t160, t163 * (-t165 + t207) - t205, 0, -t103 * t176 + t125 * t182 + t160 * t83, -t104 * t176 + t128 * t182 + t163 * t83, -t182 * t200 * t166 + t176 * t131 - t47, -t176 * t47 + t182 * t83, t163 * (-t117 * t196 + t154 * t97) - t160 * (t117 * t195 + t152 * t97), t163 * (-t152 * t81 - t154 * t79) - t160 * (-t152 * t79 + t154 * t81), t163 * (-t152 * t106 + t226) - t160 * (t154 * t106 + t227), t163 * (t115 * t195 + t152 * t177) - t160 * (t115 * t196 - t154 * t177), t163 * (t154 * t105 - t217) - t160 * (t152 * t105 + t215), (t163 * (-t115 * t154 + t117 * t152) - t160 * (-t115 * t152 - t117 * t154)) * qJD(4), t163 * (-qJ(5) * t65 + t218) - t160 * (qJ(5) * t66 - t216) + t171 * t79 - t176 * t37, t163 * (-qJ(5) * t71 + t216) - t160 * (qJ(5) * t72 + t218) + t171 * t81 - t176 * t44, t163 * (-qJ(5) * t54 - t18) - t160 * (qJ(5) * t55 + t19) + t171 * t77 - t176 * t31, t171 * t64 + t221 * t7, t163 * (t154 * t46 + t188) - t160 * (t152 * t46 - t187), t163 * (t152 * t76 + t154 * t29) - t160 * (t152 * t29 - t154 * t76), t163 * (t152 * t52 + t154 * t42) - t160 * (t152 * t42 - t154 * t52), t163 * (t154 * t45 - t188) - t160 * (t152 * t45 + t187), t163 * (-t152 * t48 + t154 * t43) - t160 * (t152 * t43 + t154 * t48), t163 * (t152 * t94 + t154 * t56) - t160 * (t152 * t56 - t154 * t94), t163 * (-qJ(5) * t22 - t152 * t9 + t154 * t16) - t160 * (qJ(5) * t23 + t152 * t16 + t154 * t9) + t171 * t38 - t176 * t11, t163 * (-qJ(5) * t24 - t152 * t10 + t154 * t17) - t160 * (qJ(5) * t25 + t154 * t10 + t152 * t17) + t171 * t40 - t176 * t12, t163 * (-qJ(5) * t20 + t154 * t4) - t160 * (qJ(5) * t21 + t152 * t4) - t176 * t8 + (pkin(5) * t203 - t160 * t185 + t182) * t28, (t163 * (pkin(5) * t152 - pkin(8) * t154) - t160 * (-pkin(8) * t152 + t185) + t182) * t5 + t221 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, 0, 0, 0, 0, 0, 0, -t160 * t135 + t163 * t137, -t163 * t134 - t160 * t138, 0, t160 * t74 + t163 * t75, 0, 0, 0, 0, 0, 0, -t160 * t65 + t163 * t66, -t160 * t71 + t163 * t72, -t160 * t54 + t163 * t55, -t160 * t18 + t163 * t19, 0, 0, 0, 0, 0, 0, -t160 * t22 + t163 * t23, -t160 * t24 + t163 * t25, -t160 * t20 + t163 * t21, -t160 * t2 + t163 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t166, t88, 0, 0, 0, 0, 0, 0, t103, t104, -t131, t47, 0, 0, 0, 0, 0, 0, t37, t44, t31, t7, 0, 0, 0, 0, 0, 0, t11, t12, t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, (-t147 + t148) * t166, t141, -t186, -t189, qJDD(4), -t74, -t75, 0, 0, t208, t114 - t113, t82, -t208, t80, qJDD(4), pkin(4) * t65 - t33, pkin(4) * t71 - t34, pkin(4) * t54, pkin(4) * t18, t102 * t209 + t159 * t70, t159 * t49 + t162 * t51, t162 * t85 + t224, t100 * t210 + t162 * t69, t159 * t84 + t211, (-t100 * t159 - t102 * t162) * t112, pkin(4) * t22 + pkin(5) * t49 + pkin(8) * t39 - t212, pkin(4) * t24 + pkin(5) * t53 + pkin(8) * t41 + t214, pkin(4) * t20 + pkin(5) * t67 + pkin(8) * t30 + t6, pkin(4) * t2 - pkin(5) * t26 + pkin(8) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t81, t77, t64, 0, 0, 0, 0, 0, 0, t38, t40, t28, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t76, t52, -t78, -t48, t94, -t14, -t15, 0, 0;];
tauJ_reg  = t13;
