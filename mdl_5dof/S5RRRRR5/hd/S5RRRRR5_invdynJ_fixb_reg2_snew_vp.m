% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:23
% EndTime: 2020-01-03 12:13:27
% DurationCPUTime: 1.38s
% Computational Cost: add. (10288->208), mult. (10766->289), div. (0->0), fcn. (6624->10), ass. (0->151)
t149 = sin(qJ(5));
t143 = qJDD(4) + qJDD(5);
t146 = qJD(1) + qJD(2);
t140 = qJD(3) + t146;
t154 = cos(qJ(5));
t155 = cos(qJ(4));
t150 = sin(qJ(4));
t195 = t140 * t150;
t106 = -t154 * t155 * t140 + t149 * t195;
t108 = (t155 * t149 + t150 * t154) * t140;
t84 = t108 * t106;
t207 = -t84 + t143;
t210 = t149 * t207;
t209 = t154 * t207;
t138 = t140 ^ 2;
t144 = qJDD(1) + qJDD(2);
t139 = qJDD(3) + t144;
t151 = sin(qJ(3));
t156 = cos(qJ(3));
t153 = sin(qJ(1));
t158 = cos(qJ(1));
t170 = -t158 * g(2) - t153 * g(3);
t127 = qJDD(1) * pkin(1) + t170;
t169 = t153 * g(2) - t158 * g(3);
t128 = -qJD(1) ^ 2 * pkin(1) - t169;
t152 = sin(qJ(2));
t157 = cos(qJ(2));
t93 = t157 * t127 - t152 * t128;
t163 = t144 * pkin(2) + t93;
t142 = t146 ^ 2;
t94 = t152 * t127 + t157 * t128;
t88 = -t142 * pkin(2) + t94;
t69 = t151 * t163 + t156 * t88;
t65 = -t138 * pkin(3) + t139 * pkin(8) + t69;
t200 = t150 * t65;
t51 = t155 * g(1) + t200;
t52 = -t150 * g(1) + t155 * t65;
t33 = t150 * t51 + t155 * t52;
t145 = qJD(4) + qJD(5);
t103 = t145 * t106;
t189 = qJD(4) * t140;
t184 = t155 * t189;
t191 = t150 * t139;
t117 = t184 + t191;
t133 = t155 * t139;
t185 = t150 * t189;
t168 = t133 - t185;
t72 = -t106 * qJD(5) + t154 * t117 + t149 * t168;
t208 = -t103 + t72;
t196 = t138 * t150;
t205 = t117 * pkin(9);
t206 = (pkin(4) * t196 + pkin(9) * t189 - g(1)) * t155 + qJDD(4) * pkin(4) - t200 - t205;
t104 = t106 ^ 2;
t105 = t108 ^ 2;
t141 = t145 ^ 2;
t126 = qJD(4) * pkin(4) - pkin(9) * t195;
t148 = t155 ^ 2;
t135 = t148 * t138;
t43 = -pkin(4) * t135 + t168 * pkin(9) - qJD(4) * t126 + t52;
t22 = t149 * t43 - t154 * t206;
t199 = t154 * t43;
t23 = t206 * t149 + t199;
t8 = t149 * t23 - t154 * t22;
t204 = t150 * t8;
t182 = t151 * t88 - t156 * t163;
t64 = -t139 * pkin(3) - t138 * pkin(8) + t182;
t203 = -pkin(3) * t64 + pkin(8) * t33;
t45 = -t168 * pkin(4) - pkin(9) * t135 + t126 * t195 + t64;
t202 = t149 * t45;
t80 = t84 + t143;
t201 = t149 * t80;
t198 = t154 * t45;
t197 = t154 * t80;
t194 = t145 * t149;
t193 = t145 * t154;
t130 = t155 * t196;
t124 = qJDD(4) + t130;
t192 = t150 * t124;
t190 = t155 * (qJDD(4) - t130);
t20 = t151 * t33 - t156 * t64;
t188 = pkin(2) * t20 + t203;
t116 = 0.2e1 * t184 + t191;
t147 = t150 ^ 2;
t134 = t147 * t138;
t159 = qJD(4) ^ 2;
t98 = -t190 - t150 * (-t134 - t159);
t187 = -pkin(3) * t116 + pkin(8) * t98 + t150 * t64;
t118 = t133 - 0.2e1 * t185;
t97 = t155 * (-t135 - t159) - t192;
t186 = pkin(3) * t118 + pkin(8) * t97 - t155 * t64;
t177 = t149 * t117 - t154 * t168;
t162 = (-qJD(5) + t145) * t108 - t177;
t62 = t103 + t72;
t34 = t149 * t162 - t154 * t62;
t35 = t149 * t62 + t154 * t162;
t14 = -t150 * t34 + t155 * t35;
t73 = -t104 - t105;
t9 = t149 * t22 + t154 * t23;
t183 = t150 * (-pkin(9) * t34 - t8) + t155 * (-pkin(4) * t73 + pkin(9) * t35 + t9) - pkin(3) * t73 + pkin(8) * t14;
t78 = -t141 - t104;
t47 = t149 * t78 + t209;
t48 = t154 * t78 - t210;
t29 = -t150 * t47 + t155 * t48;
t57 = (qJD(5) + t145) * t108 + t177;
t181 = t150 * (-pkin(9) * t47 + t202) + t155 * (-pkin(4) * t57 + pkin(9) * t48 - t198) - pkin(3) * t57 + pkin(8) * t29;
t99 = -t105 - t141;
t66 = t154 * t99 - t201;
t67 = -t149 * t99 - t197;
t37 = -t150 * t66 + t155 * t67;
t180 = t150 * (-pkin(9) * t66 + t198) + t155 * (-pkin(4) * t208 + pkin(9) * t67 + t202) - pkin(3) * t208 + pkin(8) * t37;
t77 = -t156 * t116 + t151 * t98;
t179 = pkin(2) * t77 + t187;
t76 = t156 * t118 + t151 * t97;
t178 = pkin(2) * t76 + t186;
t120 = (t147 + t148) * t139;
t123 = t134 + t135;
t176 = pkin(3) * t123 + pkin(8) * t120 + t33;
t11 = t151 * t14 - t156 * t73;
t175 = pkin(2) * t11 + t183;
t18 = t151 * t29 - t156 * t57;
t174 = pkin(2) * t18 + t181;
t25 = t151 * t37 - t156 * t208;
t173 = pkin(2) * t25 + t180;
t167 = t151 * t138 - t156 * t139;
t172 = -pkin(2) * t167 - t182;
t86 = t151 * t120 + t156 * t123;
t171 = pkin(2) * t86 + t176;
t121 = -t156 * t138 - t151 * t139;
t4 = t155 * t9 - t204;
t166 = pkin(8) * t4 - pkin(9) * t204 - pkin(3) * t45 + t155 * (-pkin(4) * t45 + pkin(9) * t9);
t2 = t151 * t4 - t156 * t45;
t165 = pkin(2) * t2 + t166;
t160 = pkin(2) * t121 - t69;
t101 = -t105 + t141;
t100 = t104 - t141;
t96 = t192 + t155 * (-t134 + t159);
t95 = t150 * (t135 - t159) + t190;
t90 = (t117 + t184) * t150;
t89 = t118 * t155;
t83 = t155 * t116 + t150 * t118;
t82 = t105 - t104;
t71 = -t108 * qJD(5) - t177;
t46 = (t150 * (-t106 * t154 + t108 * t149) + t155 * (-t106 * t149 - t108 * t154)) * t145;
t41 = t151 * t69 - t156 * t182;
t40 = pkin(2) * t41;
t39 = t150 * (t154 * t100 - t201) + t155 * (t149 * t100 + t197);
t38 = t150 * (-t149 * t101 + t209) + t155 * (t154 * t101 + t210);
t32 = t150 * (-t108 * t194 + t154 * t72) + t155 * (t108 * t193 + t149 * t72);
t31 = t150 * (t106 * t193 - t149 * t71) + t155 * (t106 * t194 + t154 * t71);
t13 = t150 * (-t149 * t208 - t154 * t57) + t155 * (-t149 * t57 + t154 * t208);
t1 = [0, 0, 0, 0, 0, qJDD(1), t170, t169, 0, 0, 0, 0, 0, 0, 0, t144, pkin(1) * (-t152 * t142 + t157 * t144) + t93, pkin(1) * (-t157 * t142 - t152 * t144) - t94, 0, pkin(1) * (t152 * t94 + t157 * t93), 0, 0, 0, 0, 0, t139, pkin(1) * (t121 * t152 - t157 * t167) + t172, pkin(1) * (t157 * t121 + t152 * t167) + t160, 0, pkin(1) * (t152 * (t151 * t182 + t156 * t69) + t157 * t41) + t40, t90, t83, t96, t89, t95, 0, pkin(1) * (t152 * (-t151 * t118 + t156 * t97) + t157 * t76) + t178, pkin(1) * (t152 * (t151 * t116 + t156 * t98) + t157 * t77) + t179, pkin(1) * (t152 * (t156 * t120 - t151 * t123) + t157 * t86) + t171, pkin(1) * (t152 * (t151 * t64 + t156 * t33) + t157 * t20) + t188, t32, t13, t38, t31, t39, t46, pkin(1) * (t152 * (t151 * t57 + t156 * t29) + t157 * t18) + t174, pkin(1) * (t152 * (t151 * t208 + t156 * t37) + t157 * t25) + t173, pkin(1) * (t152 * (t156 * t14 + t151 * t73) + t157 * t11) + t175, pkin(1) * (t152 * (t151 * t45 + t156 * t4) + t157 * t2) + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t93, -t94, 0, 0, 0, 0, 0, 0, 0, t139, t172, t160, 0, t40, t90, t83, t96, t89, t95, 0, t178, t179, t171, t188, t32, t13, t38, t31, t39, t46, t174, t173, t175, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t182, -t69, 0, 0, t90, t83, t96, t89, t95, 0, t186, t187, t176, t203, t32, t13, t38, t31, t39, t46, t181, t180, t183, t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t134 - t135, t191, t130, t133, qJDD(4), -t51, -t52, 0, 0, t84, t82, t62, -t84, t162, t143, pkin(4) * t47 - t22, -t199 - t149 * (pkin(9) * t184 - t205 - t51) + (-t149 * t124 + t66) * pkin(4), pkin(4) * t34, pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t82, t62, -t84, t162, t143, -t22, -t23, 0, 0;];
tauJ_reg = t1;
