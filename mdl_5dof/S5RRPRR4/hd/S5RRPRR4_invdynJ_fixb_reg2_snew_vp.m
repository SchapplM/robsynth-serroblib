% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:34
% EndTime: 2022-01-20 10:48:39
% DurationCPUTime: 1.44s
% Computational Cost: add. (6873->208), mult. (9072->297), div. (0->0), fcn. (5590->10), ass. (0->142)
t155 = sin(qJ(5));
t146 = qJDD(4) + qJDD(5);
t149 = qJD(1) + qJD(2);
t159 = cos(qJ(5));
t160 = cos(qJ(4));
t156 = sin(qJ(4));
t191 = t149 * t156;
t108 = -t159 * t160 * t149 + t155 * t191;
t110 = (t160 * t155 + t156 * t159) * t149;
t86 = t110 * t108;
t201 = -t86 + t146;
t204 = t155 * t201;
t203 = t159 * t201;
t152 = -g(3) + qJDD(3);
t145 = t149 ^ 2;
t147 = qJDD(1) + qJDD(2);
t153 = sin(pkin(9));
t154 = cos(pkin(9));
t158 = sin(qJ(1));
t162 = cos(qJ(1));
t181 = t158 * g(1) - t162 * g(2);
t126 = qJDD(1) * pkin(1) + t181;
t171 = t162 * g(1) + t158 * g(2);
t127 = -qJD(1) ^ 2 * pkin(1) - t171;
t157 = sin(qJ(2));
t161 = cos(qJ(2));
t91 = t161 * t126 - t157 * t127;
t166 = t147 * pkin(2) + t91;
t92 = t157 * t126 + t161 * t127;
t90 = -t145 * pkin(2) + t92;
t67 = t153 * t166 + t154 * t90;
t57 = -t145 * pkin(3) + t147 * pkin(7) + t67;
t49 = -t160 * t152 + t156 * t57;
t50 = t156 * t152 + t160 * t57;
t29 = t156 * t49 + t160 * t50;
t148 = qJD(4) + qJD(5);
t105 = t148 * t108;
t187 = qJD(4) * t149;
t182 = t160 * t187;
t189 = t156 * t147;
t119 = t182 + t189;
t137 = t160 * t147;
t183 = t156 * t187;
t170 = t137 - t183;
t72 = -t108 * qJD(5) + t159 * t119 + t155 * t170;
t202 = -t105 + t72;
t200 = -t49 + (-t119 + t182) * pkin(8);
t106 = t108 ^ 2;
t107 = t110 ^ 2;
t144 = t148 ^ 2;
t134 = t160 * t145 * t156;
t186 = qJDD(4) + t134;
t164 = t186 * pkin(4) + t200;
t130 = qJD(4) * pkin(4) - pkin(8) * t191;
t151 = t160 ^ 2;
t139 = t151 * t145;
t43 = -pkin(4) * t139 + t170 * pkin(8) - qJD(4) * t130 + t50;
t22 = t155 * t43 - t159 * t164;
t196 = t159 * t43;
t23 = t155 * t164 + t196;
t8 = t155 * t23 - t159 * t22;
t199 = t156 * t8;
t180 = t153 * t90 - t154 * t166;
t56 = -t147 * pkin(3) - t145 * pkin(7) + t180;
t45 = -t170 * pkin(4) - pkin(8) * t139 + t130 * t191 + t56;
t198 = t155 * t45;
t80 = t86 + t146;
t197 = t155 * t80;
t195 = t159 * t45;
t194 = t159 * t80;
t193 = t148 * t155;
t192 = t148 * t159;
t190 = t156 * t186;
t129 = qJDD(4) - t134;
t188 = t160 * t129;
t17 = t153 * t29 - t154 * t56;
t185 = pkin(2) * t17 - pkin(3) * t56 + pkin(7) * t29;
t122 = -t154 * t145 - t153 * t147;
t184 = pkin(2) * t122 - t67;
t9 = t155 * t22 + t159 * t23;
t150 = t156 ^ 2;
t138 = t150 * t145;
t163 = qJD(4) ^ 2;
t131 = -t138 - t163;
t100 = -t156 * t131 - t188;
t118 = 0.2e1 * t182 + t189;
t77 = t153 * t100 - t154 * t118;
t179 = pkin(2) * t77 - pkin(3) * t118 + pkin(7) * t100 + t156 * t56;
t120 = t137 - 0.2e1 * t183;
t132 = -t139 - t163;
t99 = t160 * t132 - t190;
t76 = t154 * t120 + t153 * t99;
t178 = pkin(2) * t76 + pkin(3) * t120 + pkin(7) * t99 - t160 * t56;
t177 = t155 * t119 - t159 * t170;
t165 = (-qJD(5) + t148) * t110 - t177;
t65 = t105 + t72;
t34 = t155 * t165 - t159 * t65;
t35 = t155 * t65 + t159 * t165;
t14 = -t156 * t34 + t160 * t35;
t73 = -t106 - t107;
t11 = t153 * t14 - t154 * t73;
t176 = t156 * (-pkin(8) * t34 - t8) + t160 * (-pkin(4) * t73 + pkin(8) * t35 + t9) - pkin(3) * t73 + pkin(7) * t14 + pkin(2) * t11;
t78 = -t144 - t106;
t51 = t155 * t78 + t203;
t52 = t159 * t78 - t204;
t31 = -t156 * t51 + t160 * t52;
t60 = (qJD(5) + t148) * t110 + t177;
t21 = t153 * t31 - t154 * t60;
t175 = t156 * (-pkin(8) * t51 + t198) + t160 * (-pkin(4) * t60 + pkin(8) * t52 - t195) - pkin(3) * t60 + pkin(7) * t31 + pkin(2) * t21;
t101 = -t107 - t144;
t68 = t159 * t101 - t197;
t69 = -t155 * t101 - t194;
t39 = -t156 * t68 + t160 * t69;
t25 = t153 * t39 - t154 * t202;
t174 = t156 * (-pkin(8) * t68 + t195) + t160 * (-pkin(4) * t202 + pkin(8) * t69 + t198) - pkin(3) * t202 + pkin(7) * t39 + pkin(2) * t25;
t168 = t153 * t145 - t154 * t147;
t173 = -pkin(2) * t168 - t180;
t124 = (t150 + t151) * t147;
t125 = t138 + t139;
t89 = t153 * t124 + t154 * t125;
t172 = pkin(2) * t89 + pkin(3) * t125 + pkin(7) * t124 + t29;
t4 = t160 * t9 - t199;
t2 = t153 * t4 - t154 * t45;
t167 = pkin(2) * t2 + pkin(7) * t4 - pkin(8) * t199 - pkin(3) * t45 + t160 * (-pkin(4) * t45 + pkin(8) * t9);
t103 = -t107 + t144;
t102 = t106 - t144;
t98 = t190 + t160 * (-t138 + t163);
t97 = t156 * (t139 - t163) + t188;
t94 = (t119 + t182) * t156;
t93 = t120 * t160;
t87 = t160 * t118 + t156 * t120;
t82 = t107 - t106;
t71 = -t110 * qJD(5) - t177;
t46 = (t156 * (-t108 * t159 + t110 * t155) + t160 * (-t108 * t155 - t110 * t159)) * t148;
t41 = t156 * (t159 * t102 - t197) + t160 * (t155 * t102 + t194);
t40 = t156 * (-t155 * t103 + t203) + t160 * (t159 * t103 + t204);
t37 = t153 * t67 - t154 * t180;
t36 = pkin(2) * t37;
t33 = t156 * (-t110 * t193 + t159 * t72) + t160 * (t110 * t192 + t155 * t72);
t32 = t156 * (t108 * t192 - t155 * t71) + t160 * (t108 * t193 + t159 * t71);
t13 = t156 * (-t155 * t202 - t159 * t60) + t160 * (-t155 * t60 + t159 * t202);
t1 = [0, 0, 0, 0, 0, qJDD(1), t181, t171, 0, 0, 0, 0, 0, 0, 0, t147, pkin(1) * (-t157 * t145 + t161 * t147) + t91, pkin(1) * (-t161 * t145 - t157 * t147) - t92, 0, pkin(1) * (t157 * t92 + t161 * t91), 0, 0, 0, 0, 0, t147, pkin(1) * (t157 * t122 - t161 * t168) + t173, pkin(1) * (t161 * t122 + t157 * t168) + t184, 0, pkin(1) * (t157 * (t153 * t180 + t154 * t67) + t161 * t37) + t36, t94, t87, t98, t93, t97, 0, pkin(1) * (t157 * (-t153 * t120 + t154 * t99) + t161 * t76) + t178, pkin(1) * (t157 * (t154 * t100 + t153 * t118) + t161 * t77) + t179, pkin(1) * (t157 * (t154 * t124 - t153 * t125) + t161 * t89) + t172, pkin(1) * (t157 * (t153 * t56 + t154 * t29) + t161 * t17) + t185, t33, t13, t40, t32, t41, t46, pkin(1) * (t157 * (t153 * t60 + t154 * t31) + t161 * t21) + t175, pkin(1) * (t157 * (t153 * t202 + t154 * t39) + t161 * t25) + t174, pkin(1) * (t157 * (t154 * t14 + t153 * t73) + t161 * t11) + t176, pkin(1) * (t157 * (t153 * t45 + t154 * t4) + t161 * t2) + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t91, -t92, 0, 0, 0, 0, 0, 0, 0, t147, t173, t184, 0, t36, t94, t87, t98, t93, t97, 0, t178, t179, t172, t185, t33, t13, t40, t32, t41, t46, t175, t174, t176, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, 0, 0, 0, 0, 0, 0, t156 * t132 + t160 * t186, -t156 * t129 + t160 * t131, 0, t156 * t50 - t160 * t49, 0, 0, 0, 0, 0, 0, t156 * t52 + t160 * t51, t156 * t69 + t160 * t68, t156 * t35 + t160 * t34, t156 * t9 + t160 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, t138 - t139, t189, t134, t137, qJDD(4), -t49, -t50, 0, 0, t86, t82, t65, -t86, t165, t146, pkin(4) * t51 - t22, -t196 - t155 * t200 + (-t155 * t186 + t68) * pkin(4), pkin(4) * t34, pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t82, t65, -t86, t165, t146, -t22, -t23, 0, 0;];
tauJ_reg = t1;
