% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:31:50
% EndTime: 2019-05-05 15:31:55
% DurationCPUTime: 2.42s
% Computational Cost: add. (11631->291), mult. (22106->406), div. (0->0), fcn. (14078->10), ass. (0->195)
t167 = sin(qJ(6));
t172 = cos(qJ(4));
t200 = qJD(1) * qJD(4);
t151 = t172 * t200;
t169 = sin(qJ(4));
t152 = t169 * qJDD(1);
t136 = -t152 - t151;
t130 = qJDD(5) - t136;
t127 = qJDD(6) + t130;
t168 = sin(qJ(5));
t171 = cos(qJ(5));
t204 = qJD(1) * t172;
t131 = -t171 * qJD(4) + t168 * t204;
t133 = qJD(4) * t168 + t171 * t204;
t170 = cos(qJ(6));
t111 = t170 * t131 + t133 * t167;
t113 = -t131 * t167 + t133 * t170;
t83 = t113 * t111;
t232 = -t83 + t127;
t237 = t167 * t232;
t116 = t133 * t131;
t231 = -t116 + t130;
t236 = t168 * t231;
t235 = t170 * t232;
t234 = t171 * t231;
t154 = t172 * qJDD(1);
t196 = t169 * t200;
t137 = t154 - t196;
t185 = -t168 * qJDD(4) - t171 * t137;
t107 = -qJD(5) * t131 - t185;
t193 = -t171 * qJDD(4) + t168 * t137;
t183 = qJD(5) * t133 + t193;
t64 = -t111 * qJD(6) + t170 * t107 - t167 * t183;
t150 = qJD(1) * t169 + qJD(5);
t146 = qJD(6) + t150;
t98 = t146 * t111;
t233 = -t98 + t64;
t122 = t150 * t131;
t88 = t107 + t122;
t194 = t167 * t107 + t170 * t183;
t47 = (qJD(6) - t146) * t113 + t194;
t84 = (qJD(5) - t150) * t133 + t193;
t109 = t111 ^ 2;
t110 = t113 ^ 2;
t230 = t131 ^ 2;
t129 = t133 ^ 2;
t145 = t146 ^ 2;
t149 = t150 ^ 2;
t229 = qJD(4) ^ 2;
t228 = -pkin(2) - pkin(7);
t227 = cos(qJ(1));
t226 = sin(qJ(1));
t225 = pkin(5) * t169;
t187 = -t137 + t196;
t188 = -t136 + t151;
t173 = qJD(1) ^ 2;
t156 = 0.2e1 * qJD(3) * qJD(1);
t157 = qJDD(1) * qJ(3);
t164 = sin(pkin(10));
t165 = cos(pkin(10));
t181 = t226 * g(1) - t227 * g(2);
t177 = qJDD(1) * pkin(1) + t181;
t182 = t227 * g(1) + t226 * g(2);
t178 = -t173 * pkin(1) - t182;
t206 = t164 * t177 + t165 * t178;
t191 = t156 + t157 + t206;
t94 = t228 * t173 + t191;
t70 = t188 * pkin(4) + t187 * pkin(8) + t94;
t189 = pkin(4) * t169 - pkin(8) * t172;
t134 = t189 * qJD(1);
t163 = qJDD(1) * pkin(2);
t176 = -t164 * t178 + t165 * t177;
t100 = -t173 * qJ(3) + qJDD(3) - t163 - t176;
t174 = -qJDD(1) * pkin(7) + t100;
t161 = -g(3) + qJDD(2);
t207 = t172 * t161;
t78 = -t229 * pkin(4) + qJDD(4) * pkin(8) + t207 + (-qJD(1) * t134 + t174) * t169;
t36 = t168 * t78 - t171 * t70;
t33 = t231 * pkin(5) - t88 * pkin(9) - t36;
t119 = pkin(5) * t150 - pkin(9) * t133;
t37 = t168 * t70 + t171 * t78;
t34 = -t230 * pkin(5) - t183 * pkin(9) - t150 * t119 + t37;
t13 = t167 * t34 - t170 * t33;
t14 = t167 * t33 + t170 * t34;
t7 = -t13 * t170 + t14 * t167;
t224 = t168 * t7;
t223 = t171 * t7;
t90 = t169 * t161 - t172 * t174;
t77 = -qJDD(4) * pkin(4) - t229 * pkin(8) + t134 * t204 + t90;
t39 = t183 * pkin(5) - t230 * pkin(9) + t133 * t119 + t77;
t222 = t167 * t39;
t72 = t83 + t127;
t221 = t167 * t72;
t220 = t170 * t39;
t219 = t170 * t72;
t218 = t172 * t77;
t103 = t116 + t130;
t217 = t103 * t168;
t216 = t103 * t171;
t215 = t146 * t167;
t214 = t146 * t170;
t213 = t150 * t168;
t212 = t150 * t171;
t159 = t169 ^ 2;
t211 = t159 * t173;
t160 = t172 ^ 2;
t210 = t160 * t173;
t197 = t169 * t173 * t172;
t143 = qJDD(4) + t197;
t209 = t169 * t143;
t144 = qJDD(4) - t197;
t208 = t172 * t144;
t205 = t159 + t160;
t202 = qJD(5) + t150;
t199 = t169 * t83;
t198 = t169 * t116;
t195 = pkin(1) * t164 + qJ(3);
t8 = t13 * t167 + t170 * t14;
t21 = t168 * t36 + t171 * t37;
t192 = pkin(1) * t165 - t228;
t190 = -pkin(1) * (qJDD(1) * t164 + t165 * t173) - t206;
t186 = t168 * t37 - t171 * t36;
t91 = t169 * t174 + t207;
t62 = t169 * t91 - t172 * t90;
t79 = -t145 - t109;
t40 = t167 * t79 + t235;
t184 = pkin(5) * t40 - t13;
t92 = -t110 - t145;
t53 = t170 * t92 - t221;
t180 = pkin(5) * t53 - t14;
t179 = t189 + t195;
t175 = -pkin(1) * (-qJDD(1) * t165 + t164 * t173) + t176;
t148 = -t210 - t229;
t147 = -t211 - t229;
t141 = t205 * qJDD(1);
t138 = t154 - 0.2e1 * t196;
t135 = t152 + 0.2e1 * t151;
t121 = -t129 + t149;
t120 = -t149 + t230;
t118 = t148 * t172 - t209;
t117 = t147 * t169 + t208;
t115 = t129 - t230;
t114 = -t129 - t149;
t108 = -t149 - t230;
t101 = t129 + t230;
t97 = -t173 * pkin(2) + t191;
t96 = -t110 + t145;
t95 = t109 - t145;
t89 = t202 * t131 + t185;
t87 = t107 - t122;
t85 = -t202 * t133 - t193;
t82 = t110 - t109;
t81 = -t114 * t168 - t216;
t75 = t108 * t171 - t236;
t67 = (-t111 * t170 + t113 * t167) * t146;
t66 = (-t111 * t167 - t113 * t170) * t146;
t65 = -t109 - t110;
t63 = -qJD(6) * t113 - t194;
t61 = t168 * t88 - t171 * t84;
t59 = t170 * t95 - t221;
t58 = -t167 * t96 + t235;
t57 = t167 * t95 + t219;
t56 = t170 * t96 + t237;
t55 = t169 * t81 + t172 * t89;
t54 = -t167 * t92 - t219;
t52 = t169 * t75 + t172 * t85;
t50 = t98 + t64;
t46 = (qJD(6) + t146) * t113 + t194;
t45 = -t113 * t215 + t170 * t64;
t44 = t113 * t214 + t167 * t64;
t43 = t111 * t214 - t167 * t63;
t42 = t111 * t215 + t170 * t63;
t41 = t170 * t79 - t237;
t38 = t101 * t172 + t169 * t61;
t31 = -t168 * t53 + t171 * t54;
t29 = -pkin(9) * t53 + t220;
t28 = t167 * t50 - t170 * t47;
t27 = -t167 * t233 - t170 * t46;
t26 = -t167 * t47 - t170 * t50;
t25 = -t167 * t46 + t170 * t233;
t24 = -pkin(9) * t40 + t222;
t23 = -t168 * t40 + t171 * t41;
t19 = t169 * t21 - t218;
t18 = t169 * t31 - t172 * t233;
t17 = -pkin(5) * t233 + pkin(9) * t54 + t222;
t16 = -pkin(5) * t46 + pkin(9) * t41 - t220;
t15 = t169 * t23 - t172 * t46;
t11 = -t168 * t26 + t171 * t28;
t9 = t11 * t169 - t172 * t65;
t6 = -pkin(5) * t39 + pkin(9) * t8;
t5 = -pkin(9) * t26 - t7;
t4 = -pkin(5) * t65 + pkin(9) * t28 + t8;
t3 = t171 * t8 - t224;
t1 = t169 * t3 - t172 * t39;
t2 = [0, 0, 0, 0, 0, qJDD(1), t181, t182, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t175, t190, 0, pkin(1) * (t164 * t206 + t165 * t176), qJDD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t163 + qJDD(3) - t175, t156 + 0.2e1 * t157 - t190, pkin(1) * (-t100 * t165 + t164 * t97) - pkin(2) * t100 + qJ(3) * t97, -t187 * t172, -t135 * t172 - t138 * t169, t208 - t169 * (-t210 + t229), t188 * t169, t172 * (t211 - t229) - t209, 0, -t192 * t117 + t195 * t135 + t169 * t94, -t192 * t118 + t195 * t138 + t172 * t94, -t195 * t205 * t173 + t192 * t141 - t62, -t192 * t62 + t195 * t94, t172 * (t107 * t171 - t133 * t213) + t198, t172 * (-t168 * t87 + t171 * t85) + t169 * t115, t172 * (-t121 * t168 + t234) + t169 * t88, t172 * (t131 * t212 + t168 * t183) - t198, t172 * (t120 * t171 - t217) - t169 * t84, t169 * t130 + t172 * (-t131 * t171 + t133 * t168) * t150, t168 * t218 - t169 * t36 + t179 * (t108 * t168 + t234) - t192 * t52, t171 * t218 - t169 * t37 + t179 * (t114 * t171 - t217) - t192 * t55, -t172 * t186 + t179 * (-t168 * t84 - t171 * t88) - t192 * t38, t179 * t186 - t19 * t192, t172 * (-t168 * t44 + t171 * t45) + t199, t172 * (-t168 * t25 + t171 * t27) + t169 * t82, t172 * (-t168 * t56 + t171 * t58) + t169 * t50, t172 * (-t168 * t42 + t171 * t43) - t199, t172 * (-t168 * t57 + t171 * t59) - t169 * t47, t172 * (-t168 * t66 + t171 * t67) + t169 * t127, t172 * (-t16 * t168 + t171 * t24) + t169 * t184 + t179 * (t168 * t41 + t171 * t40) - t192 * t15, t172 * (-t168 * t17 + t171 * t29) + t169 * t180 + t179 * (t168 * t54 + t171 * t53) - t192 * t18, t172 * (-t168 * t4 + t171 * t5) + t26 * t225 + t179 * (t168 * t28 + t171 * t26) - t192 * t9, t172 * (-pkin(9) * t223 - t168 * t6) + t7 * t225 + t179 * (t168 * t8 + t223) - t192 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, 0, 0, 0, 0, 0, 0, -t144 * t169 + t147 * t172, -t143 * t172 - t148 * t169, 0, t169 * t90 + t172 * t91, 0, 0, 0, 0, 0, 0, -t169 * t85 + t172 * t75, -t169 * t89 + t172 * t81, -t101 * t169 + t172 * t61, t169 * t77 + t172 * t21, 0, 0, 0, 0, 0, 0, t169 * t46 + t172 * t23, t169 * t233 + t172 * t31, t11 * t172 + t169 * t65, t169 * t39 + t172 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t173, t100, 0, 0, 0, 0, 0, 0, t117, t118, -t141, t62, 0, 0, 0, 0, 0, 0, t52, t55, t38, t19, 0, 0, 0, 0, 0, 0, t15, t18, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, (-t159 + t160) * t173, t154, -t197, -t152, qJDD(4), -t90, -t91, 0, 0, t107 * t168 + t133 * t212, t168 * t85 + t171 * t87, t121 * t171 + t236, t131 * t213 - t171 * t183, t120 * t168 + t216, (-t131 * t168 - t133 * t171) * t150, pkin(4) * t85 + pkin(8) * t75 - t171 * t77, pkin(4) * t89 + pkin(8) * t81 + t168 * t77, pkin(4) * t101 + pkin(8) * t61 + t21, -pkin(4) * t77 + pkin(8) * t21, t168 * t45 + t171 * t44, t168 * t27 + t171 * t25, t168 * t58 + t171 * t56, t168 * t43 + t171 * t42, t168 * t59 + t171 * t57, t168 * t67 + t171 * t66, -pkin(4) * t46 + pkin(8) * t23 + t16 * t171 + t168 * t24, -pkin(4) * t233 + pkin(8) * t31 + t168 * t29 + t17 * t171, -pkin(4) * t65 + pkin(8) * t11 + t168 * t5 + t171 * t4, -pkin(4) * t39 + pkin(8) * t3 - pkin(9) * t224 + t171 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t115, t88, -t116, -t84, t130, -t36, -t37, 0, 0, t83, t82, t50, -t83, -t47, t127, t184, t180, pkin(5) * t26, pkin(5) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t82, t50, -t83, -t47, t127, -t13, -t14, 0, 0;];
tauJ_reg  = t2;
