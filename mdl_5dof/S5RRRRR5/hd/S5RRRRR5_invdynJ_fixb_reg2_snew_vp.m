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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:02:06
% EndTime: 2022-01-20 12:02:12
% DurationCPUTime: 1.54s
% Computational Cost: add. (10288->208), mult. (10766->289), div. (0->0), fcn. (6624->10), ass. (0->151)
t150 = sin(qJ(5));
t144 = qJDD(4) + qJDD(5);
t147 = qJD(1) + qJD(2);
t140 = qJD(3) + t147;
t155 = cos(qJ(5));
t156 = cos(qJ(4));
t151 = sin(qJ(4));
t196 = t140 * t151;
t106 = -t155 * t156 * t140 + t150 * t196;
t108 = (t156 * t150 + t151 * t155) * t140;
t84 = t108 * t106;
t208 = -t84 + t144;
t211 = t150 * t208;
t210 = t155 * t208;
t138 = t140 ^ 2;
t145 = qJDD(1) + qJDD(2);
t139 = qJDD(3) + t145;
t152 = sin(qJ(3));
t157 = cos(qJ(3));
t154 = sin(qJ(1));
t159 = cos(qJ(1));
t183 = t154 * g(1) - t159 * g(2);
t127 = qJDD(1) * pkin(1) + t183;
t170 = t159 * g(1) + t154 * g(2);
t128 = -qJD(1) ^ 2 * pkin(1) - t170;
t153 = sin(qJ(2));
t158 = cos(qJ(2));
t93 = t158 * t127 - t153 * t128;
t164 = t145 * pkin(2) + t93;
t143 = t147 ^ 2;
t94 = t153 * t127 + t158 * t128;
t88 = -t143 * pkin(2) + t94;
t69 = t152 * t164 + t157 * t88;
t65 = -t138 * pkin(3) + t139 * pkin(8) + t69;
t201 = t151 * t65;
t51 = t156 * g(3) + t201;
t52 = -t151 * g(3) + t156 * t65;
t33 = t151 * t51 + t156 * t52;
t146 = qJD(4) + qJD(5);
t103 = t146 * t106;
t190 = qJD(4) * t140;
t185 = t156 * t190;
t192 = t151 * t139;
t117 = t185 + t192;
t133 = t156 * t139;
t186 = t151 * t190;
t169 = t133 - t186;
t72 = -t106 * qJD(5) + t155 * t117 + t150 * t169;
t209 = -t103 + t72;
t197 = t138 * t151;
t206 = t117 * pkin(9);
t207 = (pkin(4) * t197 + pkin(9) * t190 - g(3)) * t156 + qJDD(4) * pkin(4) - t201 - t206;
t104 = t106 ^ 2;
t105 = t108 ^ 2;
t142 = t146 ^ 2;
t126 = qJD(4) * pkin(4) - pkin(9) * t196;
t149 = t156 ^ 2;
t135 = t149 * t138;
t43 = -pkin(4) * t135 + t169 * pkin(9) - qJD(4) * t126 + t52;
t22 = t150 * t43 - t155 * t207;
t200 = t155 * t43;
t23 = t207 * t150 + t200;
t8 = t150 * t23 - t155 * t22;
t205 = t151 * t8;
t182 = t152 * t88 - t157 * t164;
t64 = -t139 * pkin(3) - t138 * pkin(8) + t182;
t204 = -pkin(3) * t64 + pkin(8) * t33;
t45 = -t169 * pkin(4) - pkin(9) * t135 + t126 * t196 + t64;
t203 = t150 * t45;
t80 = t84 + t144;
t202 = t150 * t80;
t199 = t155 * t45;
t198 = t155 * t80;
t195 = t146 * t150;
t194 = t146 * t155;
t130 = t156 * t197;
t124 = qJDD(4) + t130;
t193 = t151 * t124;
t191 = t156 * (qJDD(4) - t130);
t20 = t152 * t33 - t157 * t64;
t189 = pkin(2) * t20 + t204;
t116 = 0.2e1 * t185 + t192;
t148 = t151 ^ 2;
t134 = t148 * t138;
t160 = qJD(4) ^ 2;
t98 = -t191 - t151 * (-t134 - t160);
t188 = -pkin(3) * t116 + pkin(8) * t98 + t151 * t64;
t118 = t133 - 0.2e1 * t186;
t96 = t156 * (-t135 - t160) - t193;
t187 = pkin(3) * t118 + pkin(8) * t96 - t156 * t64;
t177 = t150 * t117 - t155 * t169;
t163 = (-qJD(5) + t146) * t108 - t177;
t62 = t103 + t72;
t34 = t150 * t163 - t155 * t62;
t35 = t150 * t62 + t155 * t163;
t14 = -t151 * t34 + t156 * t35;
t73 = -t104 - t105;
t9 = t150 * t22 + t155 * t23;
t184 = t151 * (-pkin(9) * t34 - t8) + t156 * (-pkin(4) * t73 + pkin(9) * t35 + t9) - pkin(3) * t73 + pkin(8) * t14;
t78 = -t142 - t104;
t47 = t150 * t78 + t210;
t48 = t155 * t78 - t211;
t29 = -t151 * t47 + t156 * t48;
t57 = (qJD(5) + t146) * t108 + t177;
t181 = t151 * (-pkin(9) * t47 + t203) + t156 * (-pkin(4) * t57 + pkin(9) * t48 - t199) - pkin(3) * t57 + pkin(8) * t29;
t99 = -t105 - t142;
t66 = t155 * t99 - t202;
t67 = -t150 * t99 - t198;
t37 = -t151 * t66 + t156 * t67;
t180 = t151 * (-pkin(9) * t66 + t199) + t156 * (-pkin(4) * t209 + pkin(9) * t67 + t203) - pkin(3) * t209 + pkin(8) * t37;
t77 = -t157 * t116 + t152 * t98;
t179 = pkin(2) * t77 + t188;
t76 = t157 * t118 + t152 * t96;
t178 = pkin(2) * t76 + t187;
t120 = (t148 + t149) * t139;
t123 = t134 + t135;
t176 = pkin(3) * t123 + pkin(8) * t120 + t33;
t11 = t152 * t14 - t157 * t73;
t175 = pkin(2) * t11 + t184;
t18 = t152 * t29 - t157 * t57;
t174 = pkin(2) * t18 + t181;
t25 = t152 * t37 - t157 * t209;
t173 = pkin(2) * t25 + t180;
t168 = t152 * t138 - t157 * t139;
t172 = -pkin(2) * t168 - t182;
t86 = t152 * t120 + t157 * t123;
t171 = pkin(2) * t86 + t176;
t121 = -t157 * t138 - t152 * t139;
t4 = t156 * t9 - t205;
t167 = pkin(8) * t4 - pkin(9) * t205 - pkin(3) * t45 + t156 * (-pkin(4) * t45 + pkin(9) * t9);
t2 = t152 * t4 - t157 * t45;
t166 = pkin(2) * t2 + t167;
t161 = pkin(2) * t121 - t69;
t101 = -t105 + t142;
t100 = t104 - t142;
t97 = t191 + t151 * (t135 - t160);
t95 = t156 * (-t134 + t160) + t193;
t90 = t118 * t156;
t89 = (t117 + t185) * t151;
t83 = t156 * t116 + t151 * t118;
t82 = t105 - t104;
t71 = -t108 * qJD(5) - t177;
t46 = (t156 * (-t106 * t150 - t108 * t155) + t151 * (-t106 * t155 + t108 * t150)) * t146;
t41 = t152 * t69 - t157 * t182;
t40 = pkin(2) * t41;
t39 = t156 * (t150 * t100 + t198) + t151 * (t155 * t100 - t202);
t38 = t156 * (t155 * t101 + t211) + t151 * (-t150 * t101 + t210);
t32 = t156 * (t108 * t194 + t150 * t72) + t151 * (-t108 * t195 + t155 * t72);
t31 = t156 * (t106 * t195 + t155 * t71) + t151 * (t106 * t194 - t150 * t71);
t13 = t156 * (-t150 * t57 + t155 * t209) + t151 * (-t150 * t209 - t155 * t57);
t1 = [0, 0, 0, 0, 0, qJDD(1), t183, t170, 0, 0, 0, 0, 0, 0, 0, t145, pkin(1) * (-t153 * t143 + t158 * t145) + t93, pkin(1) * (-t158 * t143 - t153 * t145) - t94, 0, pkin(1) * (t153 * t94 + t158 * t93), 0, 0, 0, 0, 0, t139, pkin(1) * (t121 * t153 - t158 * t168) + t172, pkin(1) * (t158 * t121 + t153 * t168) + t161, 0, pkin(1) * (t153 * (t152 * t182 + t157 * t69) + t158 * t41) + t40, t89, t83, t95, t90, t97, 0, pkin(1) * (t153 * (-t152 * t118 + t157 * t96) + t158 * t76) + t178, pkin(1) * (t153 * (t152 * t116 + t157 * t98) + t158 * t77) + t179, pkin(1) * (t153 * (t157 * t120 - t152 * t123) + t158 * t86) + t171, pkin(1) * (t153 * (t152 * t64 + t157 * t33) + t158 * t20) + t189, t32, t13, t38, t31, t39, t46, pkin(1) * (t153 * (t152 * t57 + t157 * t29) + t158 * t18) + t174, pkin(1) * (t153 * (t152 * t209 + t157 * t37) + t158 * t25) + t173, pkin(1) * (t153 * (t157 * t14 + t152 * t73) + t158 * t11) + t175, pkin(1) * (t153 * (t152 * t45 + t157 * t4) + t158 * t2) + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, t93, -t94, 0, 0, 0, 0, 0, 0, 0, t139, t172, t161, 0, t40, t89, t83, t95, t90, t97, 0, t178, t179, t171, t189, t32, t13, t38, t31, t39, t46, t174, t173, t175, t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t182, -t69, 0, 0, t89, t83, t95, t90, t97, 0, t187, t188, t176, t204, t32, t13, t38, t31, t39, t46, t181, t180, t184, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t134 - t135, t192, t130, t133, qJDD(4), -t51, -t52, 0, 0, t84, t82, t62, -t84, t163, t144, pkin(4) * t47 - t22, -t200 - t150 * (pkin(9) * t185 - t206 - t51) + (-t150 * t124 + t66) * pkin(4), pkin(4) * t34, pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t82, t62, -t84, t163, t144, -t22, -t23, 0, 0;];
tauJ_reg = t1;
