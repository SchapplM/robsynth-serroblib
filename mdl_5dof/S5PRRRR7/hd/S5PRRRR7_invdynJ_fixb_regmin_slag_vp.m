% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:13:03
% EndTime: 2019-12-05 17:13:10
% DurationCPUTime: 2.27s
% Computational Cost: add. (1915->265), mult. (4261->381), div. (0->0), fcn. (3275->14), ass. (0->158)
t124 = qJ(3) + qJ(4);
t119 = qJ(5) + t124;
t111 = sin(t119);
t112 = cos(t119);
t125 = sin(pkin(9));
t126 = cos(pkin(9));
t127 = sin(qJ(5));
t131 = cos(qJ(5));
t132 = cos(qJ(4));
t133 = cos(qJ(3));
t199 = qJD(2) * t133;
t181 = t132 * t199;
t128 = sin(qJ(4));
t129 = sin(qJ(3));
t201 = qJD(2) * t129;
t182 = t128 * t201;
t75 = -t181 + t182;
t77 = -t128 * t199 - t132 * t201;
t160 = t127 * t75 + t131 * t77;
t120 = qJDD(3) + qJDD(4);
t130 = sin(qJ(2));
t193 = t130 * qJD(1);
t100 = qJD(2) * pkin(6) + t193;
t172 = pkin(7) * qJD(2) + t100;
t70 = t172 * t133;
t57 = t132 * t70;
t69 = t172 * t129;
t58 = qJD(3) * pkin(3) - t69;
t159 = -t128 * t58 - t57;
t190 = qJD(2) * qJD(3);
t179 = t133 * t190;
t188 = t129 * qJDD(2);
t228 = t179 + t188;
t134 = cos(qJ(2));
t191 = qJD(1) * qJD(2);
t82 = qJDD(2) * pkin(6) + t130 * qJDD(1) + t134 * t191;
t31 = -t133 * t100 * qJD(3) + qJDD(3) * pkin(3) - t228 * pkin(7) - t129 * t82;
t180 = t129 * t190;
t187 = t133 * qJDD(2);
t198 = qJD(3) * t129;
t32 = -t100 * t198 + t133 * t82 + (-t180 + t187) * pkin(7);
t146 = t159 * qJD(4) - t128 * t32 + t132 * t31;
t121 = qJD(3) + qJD(4);
t28 = qJD(4) * t181 - t121 * t182 + t128 * t187 + t228 * t132;
t2 = t120 * pkin(4) - t28 * pkin(8) + t146;
t208 = t126 * t134;
t209 = t125 * t134;
t220 = g(3) * t130;
t197 = qJD(4) * t128;
t224 = -(qJD(4) * t58 + t32) * t132 - t128 * t31 + t70 * t197;
t86 = t128 * t133 + t132 * t129;
t229 = t121 * t86;
t29 = qJD(2) * t229 + t128 * t188 - t132 * t187;
t3 = -t29 * pkin(8) - t224;
t114 = -t133 * pkin(3) - pkin(2);
t192 = t134 * qJD(1);
t83 = t114 * qJD(2) - t192;
t44 = t75 * pkin(4) + t83;
t233 = -g(1) * (-t111 * t208 + t125 * t112) - g(2) * (-t111 * t209 - t126 * t112) - t127 * t3 + t131 * t2 + t111 * t220 + t44 * t160;
t222 = t75 * pkin(8);
t19 = -t159 - t222;
t195 = qJD(5) * t127;
t38 = t127 * t77 - t131 * t75;
t232 = -g(1) * (-t125 * t111 - t112 * t208) - g(2) * (t126 * t111 - t112 * t209) + t19 * t195 + t112 * t220 - t44 * t38;
t55 = t128 * t70;
t175 = t132 * t58 - t55;
t73 = t77 * pkin(8);
t18 = t175 + t73;
t16 = t121 * pkin(4) + t18;
t211 = t131 * t19;
t163 = -t127 * t16 - t211;
t231 = qJD(5) * t163 + t233;
t185 = -qJD(4) - qJD(5);
t116 = qJD(3) - t185;
t230 = (-t19 * t116 - t2) * t127 + t232;
t218 = t160 * t38;
t166 = g(1) * t126 + g(2) * t125;
t157 = t166 * t130;
t219 = g(3) * t134;
t148 = t157 - t219;
t6 = t160 ^ 2 - t38 ^ 2;
t194 = qJD(5) * t131;
t7 = -t127 * t29 + t131 * t28 - t75 * t194 + t77 * t195;
t4 = -t38 * t116 + t7;
t144 = qJD(5) * t160 - t127 * t28 - t131 * t29;
t5 = -t116 * t160 + t144;
t85 = t128 * t129 - t132 * t133;
t72 = t85 * t130;
t152 = t85 * t134;
t196 = qJD(4) * t132;
t223 = pkin(6) + pkin(7);
t183 = qJD(3) * t223;
t90 = t129 * t183;
t91 = t133 * t183;
t94 = t223 * t129;
t95 = t223 * t133;
t227 = -qJD(1) * t152 + t128 * t91 + t132 * t90 + t94 * t196 + t95 * t197;
t153 = t86 * t134;
t213 = -t128 * t94 + t132 * t95;
t226 = qJD(1) * t153 - t213 * qJD(4) + t128 * t90 - t132 * t91;
t155 = pkin(3) * t198 - t193;
t135 = qJD(3) ^ 2;
t169 = -t134 * qJDD(1) + t130 * t191;
t225 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t135 + t130 * (t166 + t191) - t169 - t219;
t221 = pkin(3) * t128;
t215 = t77 * t75;
t214 = -t132 * t69 - t55;
t212 = qJD(2) * pkin(2);
t115 = qJDD(5) + t120;
t207 = t127 * t115;
t206 = t128 * t131;
t205 = t131 * t115;
t204 = qJDD(1) - g(3);
t122 = t129 ^ 2;
t203 = -t133 ^ 2 + t122;
t136 = qJD(2) ^ 2;
t202 = t135 + t136;
t200 = qJD(2) * t130;
t189 = qJDD(3) * t129;
t186 = t134 * qJDD(2);
t178 = -qJD(5) * t16 - t3;
t174 = t128 * t69 - t57;
t173 = -t128 * t95 - t132 * t94;
t170 = pkin(4) * t229 + t155;
t33 = -t86 * pkin(8) + t173;
t168 = pkin(8) * t229 - qJD(5) * t33 + t227;
t34 = -t85 * pkin(8) + t213;
t45 = t121 * t85;
t167 = -t45 * pkin(8) + qJD(5) * t34 - t226;
t165 = g(1) * t125 - g(2) * t126;
t71 = t86 * t130;
t162 = t127 * t72 - t131 * t71;
t161 = -t127 * t71 - t131 * t72;
t42 = t127 * t86 + t131 * t85;
t43 = -t127 * t85 + t131 * t86;
t156 = t166 * t134;
t47 = pkin(3) * t180 + t114 * qJDD(2) + t169;
t101 = -t192 - t212;
t143 = -pkin(6) * qJDD(3) + (t101 + t192 - t212) * qJD(3);
t141 = -t101 * qJD(2) + t156 + t220 - t82;
t117 = sin(t124);
t118 = cos(t124);
t138 = -g(1) * (-t125 * t117 - t118 * t208) - g(2) * (t126 * t117 - t118 * t209) + t83 * t75 + t118 * t220 + t224;
t137 = -g(1) * (-t117 * t208 + t125 * t118) - g(2) * (-t117 * t209 - t126 * t118) + t83 * t77 + t146 + t117 * t220;
t113 = t132 * pkin(3) + pkin(4);
t62 = t85 * pkin(4) + t114;
t49 = pkin(3) * t201 - t77 * pkin(4);
t30 = -t75 ^ 2 + t77 ^ 2;
t23 = t73 + t214;
t22 = t174 + t222;
t21 = -qJD(2) * t153 + t121 * t72;
t20 = -qJD(2) * t152 - t130 * t229;
t15 = -t77 * t121 - t29;
t14 = t75 * t121 + t28;
t13 = t29 * pkin(4) + t47;
t10 = qJD(5) * t43 - t127 * t45 + t131 * t229;
t9 = -qJD(5) * t42 - t127 * t229 - t131 * t45;
t1 = [t204, 0, -t136 * t130 + t186, -qJDD(2) * t130 - t136 * t134, 0, 0, 0, 0, 0, (-0.2e1 * t180 + t187) * t134 + (-t202 * t133 - t189) * t130, (-qJDD(3) * t130 - 0.2e1 * t134 * t190) * t133 + (t202 * t130 - t186) * t129, 0, 0, 0, 0, 0, -t71 * t120 + t21 * t121 - t134 * t29 + t75 * t200, t72 * t120 - t20 * t121 - t134 * t28 - t77 * t200, 0, 0, 0, 0, 0, (-qJD(5) * t161 - t127 * t20 + t131 * t21) * t116 + t162 * t115 - t38 * t200 + t134 * t144, -(qJD(5) * t162 + t127 * t21 + t131 * t20) * t116 - t161 * t115 - t160 * t200 - t134 * t7; 0, qJDD(2), t204 * t134 + t157, -t204 * t130 + t156, t122 * qJDD(2) + 0.2e1 * t129 * t179, 0.2e1 * t129 * t187 - 0.2e1 * t203 * t190, t135 * t133 + t189, qJDD(3) * t133 - t135 * t129, 0, t143 * t129 + t225 * t133, -t225 * t129 + t143 * t133, t28 * t86 + t77 * t45, t229 * t77 - t28 * t85 - t86 * t29 + t45 * t75, t86 * t120 - t45 * t121, -t85 * t120 - t121 * t229, 0, t114 * t29 + t148 * t118 + t173 * t120 + t226 * t121 + t155 * t75 + t229 * t83 + t47 * t85, t114 * t28 - t117 * t148 - t213 * t120 + t227 * t121 - t155 * t77 - t83 * t45 + t47 * t86, -t160 * t9 + t7 * t43, t10 * t160 + t144 * t43 + t38 * t9 - t7 * t42, t43 * t115 + t9 * t116, -t10 * t116 - t42 * t115, 0, (-t127 * t34 + t131 * t33) * t115 - t62 * t144 + t13 * t42 + t44 * t10 - t170 * t38 + (t127 * t168 - t131 * t167) * t116 + t148 * t112, -(t127 * t33 + t131 * t34) * t115 + t62 * t7 + t13 * t43 + t44 * t9 - t170 * t160 + (t127 * t167 + t131 * t168) * t116 - t148 * t111; 0, 0, 0, 0, -t129 * t136 * t133, t203 * t136, t188, t187, qJDD(3), t129 * t141 - t133 * t165, t129 * t165 + t133 * t141, -t215, t30, t14, t15, t120, -t174 * t121 + (t132 * t120 - t121 * t197 - t75 * t201) * pkin(3) + t137, t214 * t121 + (-t128 * t120 - t121 * t196 + t77 * t201) * pkin(3) + t138, t218, t6, t4, t5, t115, t113 * t205 - (-t127 * t23 + t131 * t22) * t116 + t49 * t38 + (-t128 * t207 + (-t127 * t132 - t206) * t116 * qJD(4)) * pkin(3) + ((-pkin(3) * t206 - t127 * t113) * t116 + t163) * qJD(5) + t233, t49 * t160 + (-t113 * t115 - t2 + (-t185 * t221 + t22) * t116) * t127 + (-t115 * t221 + (-pkin(3) * t196 - qJD(5) * t113 + t23) * t116 + t178) * t131 + t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, t30, t14, t15, t120, -t121 * t159 + t137, t175 * t121 + t138, t218, t6, t4, t5, t115, -(-t127 * t18 - t211) * t116 + (-t116 * t195 - t38 * t77 + t205) * pkin(4) + t231, (t18 * t116 + t178) * t131 + (-t116 * t194 - t160 * t77 - t207) * pkin(4) + t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, t6, t4, t5, t115, -t116 * t163 + t231, (-t3 + (-qJD(5) + t116) * t16) * t131 + t230;];
tau_reg = t1;
