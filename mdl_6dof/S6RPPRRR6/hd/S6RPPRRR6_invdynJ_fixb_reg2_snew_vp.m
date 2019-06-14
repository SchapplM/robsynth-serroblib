% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:58:13
% EndTime: 2019-05-05 15:58:19
% DurationCPUTime: 2.21s
% Computational Cost: add. (9843->291), mult. (18844->375), div. (0->0), fcn. (11520->8), ass. (0->195)
t159 = sin(qJ(6));
t165 = cos(qJ(4));
t193 = qJD(1) * qJD(4);
t148 = t165 * t193;
t161 = sin(qJ(4));
t192 = t161 * qJDD(1);
t136 = -t148 - t192;
t130 = qJDD(5) - t136;
t127 = qJDD(6) + t130;
t160 = sin(qJ(5));
t164 = cos(qJ(5));
t197 = qJD(1) * t165;
t132 = -qJD(4) * t164 + t160 * t197;
t134 = qJD(4) * t160 + t164 * t197;
t163 = cos(qJ(6));
t108 = t132 * t163 + t134 * t159;
t110 = -t132 * t159 + t134 * t163;
t80 = t110 * t108;
t221 = t127 - t80;
t226 = t159 * t221;
t116 = t134 * t132;
t220 = -t116 + t130;
t225 = t160 * t220;
t224 = t163 * t221;
t223 = t164 * t220;
t186 = t161 * t193;
t191 = t165 * qJDD(1);
t137 = -t186 + t191;
t182 = -qJDD(4) * t164 + t160 * t137;
t101 = -qJD(5) * t134 - t182;
t177 = -t160 * qJDD(4) - t164 * t137;
t102 = -qJD(5) * t132 - t177;
t63 = -qJD(6) * t108 + t101 * t159 + t102 * t163;
t145 = qJD(1) * t161 + qJD(5);
t143 = qJD(6) + t145;
t93 = t143 * t108;
t222 = t63 - t93;
t123 = t145 * t132;
t88 = t102 + t123;
t184 = -t101 * t163 + t159 * t102;
t43 = (qJD(6) - t143) * t110 + t184;
t84 = (qJD(5) - t145) * t134 + t182;
t106 = t108 ^ 2;
t107 = t110 ^ 2;
t128 = t132 ^ 2;
t129 = t134 ^ 2;
t142 = t143 ^ 2;
t144 = t145 ^ 2;
t168 = qJD(1) ^ 2;
t219 = 2 * qJD(3);
t158 = t168 * pkin(7);
t162 = sin(qJ(1));
t166 = cos(qJ(1));
t185 = t162 * g(1) - g(2) * t166;
t178 = -qJDD(2) + t185;
t170 = t168 * qJ(2) + t178;
t215 = pkin(1) + qJ(3);
t183 = t215 * qJDD(1);
t169 = t183 + t170;
t79 = -t136 * pkin(4) - t137 * pkin(8) - t158 + (t219 + (pkin(4) * t165 + pkin(8) * t161) * qJD(4)) * qJD(1) + t169;
t155 = qJDD(1) * qJ(2);
t179 = t166 * g(1) + t162 * g(2);
t175 = 0.2e1 * qJD(2) * qJD(1) - t179;
t172 = qJDD(3) + t175;
t120 = -t168 * t215 + t155 + t172;
t115 = -qJDD(1) * pkin(7) + t120;
t104 = t165 * g(3) - t115 * t161;
t167 = qJD(4) ^ 2;
t181 = pkin(4) * t161 - pkin(8) * t165;
t174 = t168 * t181;
t83 = -t167 * pkin(4) + qJDD(4) * pkin(8) - t161 * t174 - t104;
t54 = t160 * t83 - t164 * t79;
t33 = pkin(5) * t220 - pkin(9) * t88 - t54;
t180 = pkin(5) * t145 - pkin(9) * t134;
t55 = t160 * t79 + t164 * t83;
t34 = -t128 * pkin(5) + t101 * pkin(9) - t145 * t180 + t55;
t15 = t159 * t34 - t163 * t33;
t16 = t159 * t33 + t163 * t34;
t7 = -t15 * t163 + t159 * t16;
t218 = t160 * t7;
t217 = t161 * pkin(5);
t216 = t164 * t7;
t214 = qJ(2) - pkin(7);
t103 = g(3) * t161 + t115 * t165;
t82 = qJDD(4) * pkin(4) + pkin(8) * t167 - t165 * t174 + t103;
t48 = pkin(5) * t101 + pkin(9) * t128 - t134 * t180 + t82;
t213 = t159 * t48;
t68 = t127 + t80;
t212 = t159 * t68;
t97 = t116 + t130;
t211 = t160 * t97;
t210 = t163 * t48;
t209 = t163 * t68;
t208 = t164 * t97;
t207 = t165 * t82;
t206 = qJDD(1) * pkin(1);
t205 = t143 * t159;
t204 = t143 * t163;
t203 = t145 * t160;
t202 = t145 * t164;
t156 = t161 ^ 2;
t201 = t156 * t168;
t157 = t165 ^ 2;
t200 = t157 * t168;
t187 = t161 * t168 * t165;
t199 = t161 * (qJDD(4) + t187);
t198 = t156 + t157;
t195 = qJD(5) + t145;
t190 = qJD(1) * t219;
t189 = t161 * t80;
t188 = t161 * t116;
t8 = t15 * t159 + t16 * t163;
t31 = t160 * t54 + t164 * t55;
t30 = t160 * t55 - t164 * t54;
t72 = t165 * t103 - t161 * t104;
t73 = -t142 - t106;
t36 = t159 * t73 + t224;
t176 = pkin(5) * t36 - t15;
t173 = t181 + t215;
t90 = -t107 - t142;
t50 = t163 * t90 - t212;
t171 = pkin(5) * t50 - t16;
t119 = t169 + t190;
t151 = 0.2e1 * t155;
t140 = t198 * t168;
t139 = t198 * qJDD(1);
t138 = -0.2e1 * t186 + t191;
t135 = 0.2e1 * t148 + t192;
t131 = t165 * (qJDD(4) - t187);
t124 = t170 + t206;
t122 = -t129 + t144;
t121 = t128 - t144;
t118 = -t199 + t165 * (-t167 - t200);
t117 = t161 * (-t167 - t201) + t131;
t114 = -t158 + t119;
t113 = t129 - t128;
t112 = -t129 - t144;
t105 = -t144 - t128;
t95 = t128 + t129;
t92 = -t107 + t142;
t91 = t106 - t142;
t89 = t132 * t195 + t177;
t87 = t102 - t123;
t85 = -t134 * t195 - t182;
t78 = t107 - t106;
t75 = -t112 * t160 - t208;
t74 = t112 * t164 - t211;
t71 = t105 * t164 - t225;
t70 = t105 * t160 + t223;
t66 = (-t108 * t163 + t110 * t159) * t143;
t65 = (-t108 * t159 - t110 * t163) * t143;
t64 = -t106 - t107;
t62 = -qJD(6) * t110 - t184;
t61 = t160 * t88 - t164 * t84;
t60 = -t160 * t84 - t164 * t88;
t59 = t163 * t91 - t212;
t58 = -t159 * t92 + t224;
t57 = t159 * t91 + t209;
t56 = t163 * t92 + t226;
t52 = t161 * t75 + t165 * t89;
t51 = -t159 * t90 - t209;
t49 = t161 * t71 + t165 * t85;
t46 = t63 + t93;
t42 = (qJD(6) + t143) * t110 + t184;
t41 = -t110 * t205 + t163 * t63;
t40 = t110 * t204 + t159 * t63;
t39 = t108 * t204 - t159 * t62;
t38 = t108 * t205 + t163 * t62;
t37 = t163 * t73 - t226;
t35 = t161 * t61 + t165 * t95;
t29 = -t160 * t50 + t164 * t51;
t28 = t160 * t51 + t164 * t50;
t27 = -pkin(9) * t50 - t210;
t26 = -pkin(9) * t36 - t213;
t25 = t159 * t46 - t163 * t43;
t24 = -t159 * t222 - t163 * t42;
t23 = -t159 * t43 - t163 * t46;
t22 = -t159 * t42 + t163 * t222;
t21 = -t160 * t36 + t164 * t37;
t20 = t160 * t37 + t164 * t36;
t19 = t161 * t31 + t207;
t18 = -pkin(5) * t222 + pkin(9) * t51 - t213;
t17 = t161 * t29 - t165 * t222;
t13 = -pkin(5) * t42 + pkin(9) * t37 + t210;
t12 = t161 * t21 - t165 * t42;
t11 = -t160 * t23 + t164 * t25;
t10 = t160 * t25 + t164 * t23;
t9 = t11 * t161 - t165 * t64;
t6 = pkin(5) * t48 + pkin(9) * t8;
t5 = -pkin(9) * t23 - t7;
t4 = -pkin(5) * t64 + pkin(9) * t25 + t8;
t3 = t164 * t8 - t218;
t2 = t160 * t8 + t216;
t1 = t161 * t3 + t165 * t48;
t14 = [0, 0, 0, 0, 0, qJDD(1), t185, t179, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t178 - 0.2e1 * t206, t151 + t175, qJ(2) * (-pkin(1) * t168 + t155 + t175) + pkin(1) * t124, qJDD(1), 0, 0, 0, 0, 0, 0, t151 + t172, t178 + 0.2e1 * t183 + t190, qJ(2) * t120 + t119 * t215, (t137 - t186) * t165, -t135 * t165 - t138 * t161, t131 - t161 * (t167 - t200), (-t136 + t148) * t161, t165 * (-t167 + t201) - t199, 0, t161 * t114 + t117 * t214 + t135 * t215, t165 * t114 + t118 * t214 + t138 * t215, -t139 * t214 - t140 * t215 - t72, t114 * t215 + t214 * t72, t165 * (t102 * t164 - t134 * t203) + t188, t165 * (-t160 * t87 + t164 * t85) + t161 * t113, t165 * (-t122 * t160 + t223) + t161 * t88, t165 * (-t101 * t160 + t132 * t202) - t188, t165 * (t121 * t164 - t211) - t161 * t84, t161 * t130 + t165 * (-t132 * t164 + t134 * t160) * t145, -t160 * t207 - t161 * t54 + t173 * t70 + t214 * t49, -t161 * t55 - t164 * t207 + t173 * t74 + t214 * t52, -t165 * t30 + t173 * t60 + t214 * t35, t173 * t30 + t19 * t214, t165 * (-t160 * t40 + t164 * t41) + t189, t165 * (-t160 * t22 + t164 * t24) + t161 * t78, t165 * (-t160 * t56 + t164 * t58) + t161 * t46, t165 * (-t160 * t38 + t164 * t39) - t189, t165 * (-t160 * t57 + t164 * t59) - t161 * t43, t165 * (-t160 * t65 + t164 * t66) + t161 * t127, t165 * (-t13 * t160 + t164 * t26) + t161 * t176 + t173 * t20 + t214 * t12, t165 * (-t160 * t18 + t164 * t27) + t161 * t171 + t173 * t28 + t214 * t17, t165 * (-t160 * t4 + t164 * t5) + t23 * t217 + t214 * t9 + t173 * t10, t165 * (-pkin(9) * t216 - t160 * t6) + t7 * t217 + t173 * t2 + t214 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t168, -t124, 0, 0, 0, 0, 0, 0, 0, -t168, -qJDD(1), -t119, 0, 0, 0, 0, 0, 0, -t135, -t138, t140, -t114, 0, 0, 0, 0, 0, 0, -t70, -t74, -t60, -t30, 0, 0, 0, 0, 0, 0, -t20, -t28, -t10, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t168, t120, 0, 0, 0, 0, 0, 0, t117, t118, -t139, t72, 0, 0, 0, 0, 0, 0, t49, t52, t35, t19, 0, 0, 0, 0, 0, 0, t12, t17, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, (-t156 + t157) * t168, t191, -t187, -t192, qJDD(4), t103, t104, 0, 0, t102 * t160 + t134 * t202, t160 * t85 + t164 * t87, t122 * t164 + t225, t101 * t164 + t132 * t203, t121 * t160 + t208, (-t132 * t160 - t134 * t164) * t145, pkin(4) * t85 + pkin(8) * t71 + t164 * t82, pkin(4) * t89 + pkin(8) * t75 - t160 * t82, pkin(4) * t95 + pkin(8) * t61 + t31, pkin(4) * t82 + pkin(8) * t31, t160 * t41 + t164 * t40, t160 * t24 + t164 * t22, t160 * t58 + t164 * t56, t160 * t39 + t164 * t38, t160 * t59 + t164 * t57, t160 * t66 + t164 * t65, -pkin(4) * t42 + pkin(8) * t21 + t13 * t164 + t160 * t26, -pkin(4) * t222 + pkin(8) * t29 + t160 * t27 + t164 * t18, -pkin(4) * t64 + pkin(8) * t11 + t160 * t5 + t164 * t4, pkin(4) * t48 + pkin(8) * t3 - pkin(9) * t218 + t164 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t113, t88, -t116, -t84, t130, -t54, -t55, 0, 0, t80, t78, t46, -t80, -t43, t127, t176, t171, pkin(5) * t23, pkin(5) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t78, t46, -t80, -t43, t127, -t15, -t16, 0, 0;];
tauJ_reg  = t14;
