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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:58:31
% EndTime: 2019-12-05 18:58:36
% DurationCPUTime: 1.42s
% Computational Cost: add. (10288->208), mult. (10766->289), div. (0->0), fcn. (6624->10), ass. (0->151)
t151 = sin(qJ(5));
t145 = qJDD(4) + qJDD(5);
t148 = qJD(1) + qJD(2);
t140 = qJD(3) + t148;
t156 = cos(qJ(5));
t157 = cos(qJ(4));
t152 = sin(qJ(4));
t197 = t140 * t152;
t106 = -t156 * t157 * t140 + t151 * t197;
t108 = (t157 * t151 + t152 * t156) * t140;
t84 = t108 * t106;
t209 = -t84 + t145;
t212 = t151 * t209;
t211 = t156 * t209;
t138 = t140 ^ 2;
t146 = qJDD(1) + qJDD(2);
t139 = qJDD(3) + t146;
t153 = sin(qJ(3));
t158 = cos(qJ(3));
t155 = sin(qJ(1));
t160 = cos(qJ(1));
t191 = t160 * g(2) + t155 * g(3);
t127 = qJDD(1) * pkin(1) + t191;
t171 = -t155 * g(2) + t160 * g(3);
t128 = -qJD(1) ^ 2 * pkin(1) - t171;
t154 = sin(qJ(2));
t159 = cos(qJ(2));
t93 = t159 * t127 - t154 * t128;
t165 = t146 * pkin(2) + t93;
t144 = t148 ^ 2;
t94 = t154 * t127 + t159 * t128;
t88 = -t144 * pkin(2) + t94;
t69 = t153 * t165 + t158 * t88;
t65 = -t138 * pkin(3) + t139 * pkin(8) + t69;
t202 = t152 * t65;
t51 = t157 * g(1) + t202;
t52 = -t152 * g(1) + t157 * t65;
t33 = t152 * t51 + t157 * t52;
t147 = qJD(4) + qJD(5);
t103 = t147 * t106;
t190 = qJD(4) * t140;
t185 = t157 * t190;
t193 = t152 * t139;
t117 = t185 + t193;
t133 = t157 * t139;
t186 = t152 * t190;
t170 = t133 - t186;
t72 = -t106 * qJD(5) + t156 * t117 + t151 * t170;
t210 = -t103 + t72;
t198 = t138 * t152;
t207 = t117 * pkin(9);
t208 = (pkin(4) * t198 + pkin(9) * t190 - g(1)) * t157 + qJDD(4) * pkin(4) - t202 - t207;
t104 = t106 ^ 2;
t105 = t108 ^ 2;
t143 = t147 ^ 2;
t126 = qJD(4) * pkin(4) - pkin(9) * t197;
t150 = t157 ^ 2;
t135 = t150 * t138;
t43 = -pkin(4) * t135 + pkin(9) * t170 - qJD(4) * t126 + t52;
t22 = t151 * t43 - t156 * t208;
t201 = t156 * t43;
t23 = t208 * t151 + t201;
t8 = t151 * t23 - t156 * t22;
t206 = t152 * t8;
t183 = t153 * t88 - t158 * t165;
t64 = -t139 * pkin(3) - t138 * pkin(8) + t183;
t205 = -pkin(3) * t64 + pkin(8) * t33;
t45 = -pkin(4) * t170 - pkin(9) * t135 + t126 * t197 + t64;
t204 = t151 * t45;
t80 = t84 + t145;
t203 = t151 * t80;
t200 = t156 * t45;
t199 = t156 * t80;
t196 = t147 * t151;
t195 = t147 * t156;
t130 = t157 * t198;
t124 = qJDD(4) + t130;
t194 = t152 * t124;
t192 = t157 * (qJDD(4) - t130);
t20 = t153 * t33 - t158 * t64;
t189 = pkin(2) * t20 + t205;
t116 = 0.2e1 * t185 + t193;
t149 = t152 ^ 2;
t134 = t149 * t138;
t161 = qJD(4) ^ 2;
t98 = -t192 - t152 * (-t134 - t161);
t188 = -pkin(3) * t116 + pkin(8) * t98 + t152 * t64;
t118 = t133 - 0.2e1 * t186;
t97 = t157 * (-t135 - t161) - t194;
t187 = pkin(3) * t118 + pkin(8) * t97 - t157 * t64;
t178 = t151 * t117 - t156 * t170;
t164 = (-qJD(5) + t147) * t108 - t178;
t62 = t103 + t72;
t34 = t151 * t164 - t156 * t62;
t35 = t151 * t62 + t156 * t164;
t14 = -t152 * t34 + t157 * t35;
t73 = -t104 - t105;
t9 = t151 * t22 + t156 * t23;
t184 = t152 * (-pkin(9) * t34 - t8) + t157 * (-pkin(4) * t73 + pkin(9) * t35 + t9) - pkin(3) * t73 + pkin(8) * t14;
t78 = -t143 - t104;
t47 = t151 * t78 + t211;
t48 = t156 * t78 - t212;
t29 = -t152 * t47 + t157 * t48;
t57 = (qJD(5) + t147) * t108 + t178;
t182 = t152 * (-pkin(9) * t47 + t204) + t157 * (-pkin(4) * t57 + pkin(9) * t48 - t200) - pkin(3) * t57 + pkin(8) * t29;
t99 = -t105 - t143;
t66 = t156 * t99 - t203;
t67 = -t151 * t99 - t199;
t37 = -t152 * t66 + t157 * t67;
t181 = t152 * (-pkin(9) * t66 + t200) + t157 * (-pkin(4) * t210 + pkin(9) * t67 + t204) - pkin(3) * t210 + pkin(8) * t37;
t77 = -t158 * t116 + t153 * t98;
t180 = pkin(2) * t77 + t188;
t76 = t158 * t118 + t153 * t97;
t179 = pkin(2) * t76 + t187;
t120 = (t149 + t150) * t139;
t123 = t134 + t135;
t177 = pkin(3) * t123 + pkin(8) * t120 + t33;
t11 = t153 * t14 - t158 * t73;
t176 = pkin(2) * t11 + t184;
t18 = t153 * t29 - t158 * t57;
t175 = pkin(2) * t18 + t182;
t25 = t153 * t37 - t158 * t210;
t174 = pkin(2) * t25 + t181;
t169 = t153 * t138 - t158 * t139;
t173 = -pkin(2) * t169 - t183;
t86 = t153 * t120 + t158 * t123;
t172 = pkin(2) * t86 + t177;
t121 = -t158 * t138 - t153 * t139;
t4 = t157 * t9 - t206;
t168 = pkin(8) * t4 - pkin(9) * t206 - pkin(3) * t45 + t157 * (-pkin(4) * t45 + pkin(9) * t9);
t2 = t153 * t4 - t158 * t45;
t167 = pkin(2) * t2 + t168;
t162 = pkin(2) * t121 - t69;
t101 = -t105 + t143;
t100 = t104 - t143;
t96 = t194 + t157 * (-t134 + t161);
t95 = t152 * (t135 - t161) + t192;
t90 = (t117 + t185) * t152;
t89 = t118 * t157;
t83 = t157 * t116 + t152 * t118;
t82 = t105 - t104;
t71 = -t108 * qJD(5) - t178;
t46 = (t152 * (-t106 * t156 + t108 * t151) + t157 * (-t106 * t151 - t108 * t156)) * t147;
t41 = t153 * t69 - t158 * t183;
t40 = pkin(2) * t41;
t39 = t152 * (t156 * t100 - t203) + t157 * (t151 * t100 + t199);
t38 = t152 * (-t151 * t101 + t211) + t157 * (t156 * t101 + t212);
t32 = t152 * (-t108 * t196 + t156 * t72) + t157 * (t108 * t195 + t151 * t72);
t31 = t152 * (t106 * t195 - t151 * t71) + t157 * (t106 * t196 + t156 * t71);
t13 = t152 * (-t151 * t210 - t156 * t57) + t157 * (-t151 * t57 + t156 * t210);
t1 = [0, 0, 0, 0, 0, qJDD(1), t191, t171, 0, 0, 0, 0, 0, 0, 0, t146, pkin(1) * (-t154 * t144 + t159 * t146) + t93, pkin(1) * (-t159 * t144 - t154 * t146) - t94, 0, pkin(1) * (t154 * t94 + t159 * t93), 0, 0, 0, 0, 0, t139, pkin(1) * (t121 * t154 - t159 * t169) + t173, pkin(1) * (t159 * t121 + t154 * t169) + t162, 0, pkin(1) * (t154 * (t153 * t183 + t158 * t69) + t159 * t41) + t40, t90, t83, t96, t89, t95, 0, pkin(1) * (t154 * (-t153 * t118 + t158 * t97) + t159 * t76) + t179, pkin(1) * (t154 * (t153 * t116 + t158 * t98) + t159 * t77) + t180, pkin(1) * (t154 * (t158 * t120 - t153 * t123) + t159 * t86) + t172, pkin(1) * (t154 * (t153 * t64 + t158 * t33) + t159 * t20) + t189, t32, t13, t38, t31, t39, t46, pkin(1) * (t154 * (t153 * t57 + t158 * t29) + t159 * t18) + t175, pkin(1) * (t154 * (t153 * t210 + t158 * t37) + t159 * t25) + t174, pkin(1) * (t154 * (t158 * t14 + t153 * t73) + t159 * t11) + t176, pkin(1) * (t154 * (t153 * t45 + t158 * t4) + t159 * t2) + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t93, -t94, 0, 0, 0, 0, 0, 0, 0, t139, t173, t162, 0, t40, t90, t83, t96, t89, t95, 0, t179, t180, t172, t189, t32, t13, t38, t31, t39, t46, t175, t174, t176, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t183, -t69, 0, 0, t90, t83, t96, t89, t95, 0, t187, t188, t177, t205, t32, t13, t38, t31, t39, t46, t182, t181, t184, t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t134 - t135, t193, t130, t133, qJDD(4), -t51, -t52, 0, 0, t84, t82, t62, -t84, t164, t145, pkin(4) * t47 - t22, -t201 - t151 * (pkin(9) * t185 - t207 - t51) + (-t151 * t124 + t66) * pkin(4), pkin(4) * t34, pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t82, t62, -t84, t164, t145, -t22, -t23, 0, 0;];
tauJ_reg = t1;
