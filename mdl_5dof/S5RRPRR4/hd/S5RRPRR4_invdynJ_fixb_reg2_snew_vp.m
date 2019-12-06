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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:32:09
% EndTime: 2019-12-05 18:32:13
% DurationCPUTime: 1.31s
% Computational Cost: add. (6873->208), mult. (9072->297), div. (0->0), fcn. (5590->10), ass. (0->142)
t156 = sin(qJ(5));
t147 = qJDD(4) + qJDD(5);
t150 = qJD(1) + qJD(2);
t160 = cos(qJ(5));
t161 = cos(qJ(4));
t157 = sin(qJ(4));
t192 = t150 * t157;
t108 = -t160 * t161 * t150 + t156 * t192;
t110 = (t161 * t156 + t157 * t160) * t150;
t86 = t110 * t108;
t202 = -t86 + t147;
t205 = t156 * t202;
t204 = t160 * t202;
t153 = -g(1) + qJDD(3);
t146 = t150 ^ 2;
t148 = qJDD(1) + qJDD(2);
t154 = sin(pkin(9));
t155 = cos(pkin(9));
t159 = sin(qJ(1));
t163 = cos(qJ(1));
t188 = t163 * g(2) + t159 * g(3);
t126 = qJDD(1) * pkin(1) + t188;
t172 = -t159 * g(2) + t163 * g(3);
t127 = -qJD(1) ^ 2 * pkin(1) - t172;
t158 = sin(qJ(2));
t162 = cos(qJ(2));
t91 = t162 * t126 - t158 * t127;
t167 = t148 * pkin(2) + t91;
t92 = t158 * t126 + t162 * t127;
t90 = -t146 * pkin(2) + t92;
t67 = t154 * t167 + t155 * t90;
t57 = -t146 * pkin(3) + t148 * pkin(7) + t67;
t49 = -t161 * t153 + t157 * t57;
t50 = t157 * t153 + t161 * t57;
t29 = t157 * t49 + t161 * t50;
t149 = qJD(4) + qJD(5);
t105 = t149 * t108;
t187 = qJD(4) * t150;
t182 = t161 * t187;
t190 = t157 * t148;
t119 = t182 + t190;
t137 = t161 * t148;
t183 = t157 * t187;
t171 = t137 - t183;
t72 = -t108 * qJD(5) + t160 * t119 + t156 * t171;
t203 = -t105 + t72;
t201 = -t49 + (-t119 + t182) * pkin(8);
t106 = t108 ^ 2;
t107 = t110 ^ 2;
t145 = t149 ^ 2;
t134 = t161 * t146 * t157;
t186 = qJDD(4) + t134;
t165 = t186 * pkin(4) + t201;
t130 = qJD(4) * pkin(4) - pkin(8) * t192;
t152 = t161 ^ 2;
t139 = t152 * t146;
t43 = -pkin(4) * t139 + t171 * pkin(8) - qJD(4) * t130 + t50;
t22 = t156 * t43 - t160 * t165;
t197 = t160 * t43;
t23 = t156 * t165 + t197;
t8 = t156 * t23 - t160 * t22;
t200 = t157 * t8;
t181 = t154 * t90 - t155 * t167;
t56 = -t148 * pkin(3) - t146 * pkin(7) + t181;
t45 = -t171 * pkin(4) - pkin(8) * t139 + t130 * t192 + t56;
t199 = t156 * t45;
t80 = t86 + t147;
t198 = t156 * t80;
t196 = t160 * t45;
t195 = t160 * t80;
t194 = t149 * t156;
t193 = t149 * t160;
t191 = t157 * t186;
t129 = qJDD(4) - t134;
t189 = t161 * t129;
t17 = t154 * t29 - t155 * t56;
t185 = pkin(2) * t17 - pkin(3) * t56 + pkin(7) * t29;
t122 = -t155 * t146 - t154 * t148;
t184 = pkin(2) * t122 - t67;
t9 = t156 * t22 + t160 * t23;
t151 = t157 ^ 2;
t138 = t151 * t146;
t164 = qJD(4) ^ 2;
t131 = -t138 - t164;
t100 = -t157 * t131 - t189;
t118 = 0.2e1 * t182 + t190;
t77 = t154 * t100 - t155 * t118;
t180 = pkin(2) * t77 - pkin(3) * t118 + pkin(7) * t100 + t157 * t56;
t120 = t137 - 0.2e1 * t183;
t132 = -t139 - t164;
t99 = t161 * t132 - t191;
t76 = t155 * t120 + t154 * t99;
t179 = pkin(2) * t76 + pkin(3) * t120 + pkin(7) * t99 - t161 * t56;
t178 = t156 * t119 - t160 * t171;
t166 = (-qJD(5) + t149) * t110 - t178;
t65 = t105 + t72;
t34 = t156 * t166 - t160 * t65;
t35 = t156 * t65 + t160 * t166;
t14 = -t157 * t34 + t161 * t35;
t73 = -t106 - t107;
t11 = t154 * t14 - t155 * t73;
t177 = t157 * (-pkin(8) * t34 - t8) + t161 * (-pkin(4) * t73 + pkin(8) * t35 + t9) - pkin(3) * t73 + pkin(7) * t14 + pkin(2) * t11;
t78 = -t145 - t106;
t51 = t156 * t78 + t204;
t52 = t160 * t78 - t205;
t31 = -t157 * t51 + t161 * t52;
t60 = (qJD(5) + t149) * t110 + t178;
t21 = t154 * t31 - t155 * t60;
t176 = t157 * (-pkin(8) * t51 + t199) + t161 * (-pkin(4) * t60 + pkin(8) * t52 - t196) - pkin(3) * t60 + pkin(7) * t31 + pkin(2) * t21;
t101 = -t107 - t145;
t68 = t160 * t101 - t198;
t69 = -t156 * t101 - t195;
t39 = -t157 * t68 + t161 * t69;
t25 = t154 * t39 - t155 * t203;
t175 = t157 * (-pkin(8) * t68 + t196) + t161 * (-pkin(4) * t203 + pkin(8) * t69 + t199) - pkin(3) * t203 + pkin(7) * t39 + pkin(2) * t25;
t169 = t154 * t146 - t155 * t148;
t174 = -pkin(2) * t169 - t181;
t124 = (t151 + t152) * t148;
t125 = t138 + t139;
t89 = t154 * t124 + t155 * t125;
t173 = pkin(2) * t89 + pkin(3) * t125 + pkin(7) * t124 + t29;
t4 = t161 * t9 - t200;
t2 = t154 * t4 - t155 * t45;
t168 = pkin(2) * t2 + pkin(7) * t4 - pkin(8) * t200 - pkin(3) * t45 + t161 * (-pkin(4) * t45 + pkin(8) * t9);
t103 = -t107 + t145;
t102 = t106 - t145;
t98 = t191 + t161 * (-t138 + t164);
t97 = t157 * (t139 - t164) + t189;
t94 = (t119 + t182) * t157;
t93 = t120 * t161;
t87 = t161 * t118 + t157 * t120;
t82 = t107 - t106;
t71 = -t110 * qJD(5) - t178;
t46 = (t157 * (-t108 * t160 + t110 * t156) + t161 * (-t108 * t156 - t110 * t160)) * t149;
t41 = t157 * (t160 * t102 - t198) + t161 * (t156 * t102 + t195);
t40 = t157 * (-t156 * t103 + t204) + t161 * (t160 * t103 + t205);
t37 = t154 * t67 - t155 * t181;
t36 = pkin(2) * t37;
t33 = t157 * (-t110 * t194 + t160 * t72) + t161 * (t110 * t193 + t156 * t72);
t32 = t157 * (t108 * t193 - t156 * t71) + t161 * (t108 * t194 + t160 * t71);
t13 = t157 * (-t156 * t203 - t160 * t60) + t161 * (-t156 * t60 + t160 * t203);
t1 = [0, 0, 0, 0, 0, qJDD(1), t188, t172, 0, 0, 0, 0, 0, 0, 0, t148, pkin(1) * (-t158 * t146 + t162 * t148) + t91, pkin(1) * (-t162 * t146 - t158 * t148) - t92, 0, pkin(1) * (t158 * t92 + t162 * t91), 0, 0, 0, 0, 0, t148, pkin(1) * (t158 * t122 - t162 * t169) + t174, pkin(1) * (t162 * t122 + t158 * t169) + t184, 0, pkin(1) * (t158 * (t154 * t181 + t155 * t67) + t162 * t37) + t36, t94, t87, t98, t93, t97, 0, pkin(1) * (t158 * (-t154 * t120 + t155 * t99) + t162 * t76) + t179, pkin(1) * (t158 * (t155 * t100 + t154 * t118) + t162 * t77) + t180, pkin(1) * (t158 * (t155 * t124 - t154 * t125) + t162 * t89) + t173, pkin(1) * (t158 * (t154 * t56 + t155 * t29) + t162 * t17) + t185, t33, t13, t40, t32, t41, t46, pkin(1) * (t158 * (t154 * t60 + t155 * t31) + t162 * t21) + t176, pkin(1) * (t158 * (t154 * t203 + t155 * t39) + t162 * t25) + t175, pkin(1) * (t158 * (t155 * t14 + t154 * t73) + t162 * t11) + t177, pkin(1) * (t158 * (t154 * t45 + t155 * t4) + t162 * t2) + t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t91, -t92, 0, 0, 0, 0, 0, 0, 0, t148, t174, t184, 0, t36, t94, t87, t98, t93, t97, 0, t179, t180, t173, t185, t33, t13, t40, t32, t41, t46, t176, t175, t177, t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, 0, 0, 0, 0, 0, 0, t157 * t132 + t161 * t186, -t157 * t129 + t161 * t131, 0, t157 * t50 - t161 * t49, 0, 0, 0, 0, 0, 0, t157 * t52 + t161 * t51, t157 * t69 + t161 * t68, t157 * t35 + t161 * t34, t157 * t9 + t161 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, t138 - t139, t190, t134, t137, qJDD(4), -t49, -t50, 0, 0, t86, t82, t65, -t86, t166, t147, pkin(4) * t51 - t22, -t197 - t156 * t201 + (-t156 * t186 + t68) * pkin(4), pkin(4) * t34, pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t82, t65, -t86, t166, t147, -t22, -t23, 0, 0;];
tauJ_reg = t1;
