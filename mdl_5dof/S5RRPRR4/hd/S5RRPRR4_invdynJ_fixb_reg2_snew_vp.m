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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:02:07
% EndTime: 2020-01-03 12:02:11
% DurationCPUTime: 1.30s
% Computational Cost: add. (6873->208), mult. (9072->297), div. (0->0), fcn. (5590->10), ass. (0->142)
t154 = sin(qJ(5));
t145 = qJDD(4) + qJDD(5);
t148 = qJD(1) + qJD(2);
t158 = cos(qJ(5));
t159 = cos(qJ(4));
t155 = sin(qJ(4));
t190 = t148 * t155;
t108 = -t158 * t159 * t148 + t154 * t190;
t110 = (t159 * t154 + t155 * t158) * t148;
t86 = t110 * t108;
t200 = -t86 + t145;
t203 = t154 * t200;
t202 = t158 * t200;
t151 = -g(1) + qJDD(3);
t144 = t148 ^ 2;
t146 = qJDD(1) + qJDD(2);
t152 = sin(pkin(9));
t153 = cos(pkin(9));
t157 = sin(qJ(1));
t161 = cos(qJ(1));
t171 = -t161 * g(2) - t157 * g(3);
t126 = qJDD(1) * pkin(1) + t171;
t170 = t157 * g(2) - t161 * g(3);
t127 = -qJD(1) ^ 2 * pkin(1) - t170;
t156 = sin(qJ(2));
t160 = cos(qJ(2));
t91 = t160 * t126 - t156 * t127;
t165 = t146 * pkin(2) + t91;
t92 = t156 * t126 + t160 * t127;
t90 = -t144 * pkin(2) + t92;
t67 = t152 * t165 + t153 * t90;
t57 = -t144 * pkin(3) + t146 * pkin(7) + t67;
t49 = -t159 * t151 + t155 * t57;
t50 = t155 * t151 + t159 * t57;
t29 = t155 * t49 + t159 * t50;
t147 = qJD(4) + qJD(5);
t105 = t147 * t108;
t186 = qJD(4) * t148;
t181 = t159 * t186;
t188 = t155 * t146;
t119 = t181 + t188;
t137 = t159 * t146;
t182 = t155 * t186;
t169 = t137 - t182;
t72 = -t108 * qJD(5) + t158 * t119 + t154 * t169;
t201 = -t105 + t72;
t199 = -t49 + (-t119 + t181) * pkin(8);
t106 = t108 ^ 2;
t107 = t110 ^ 2;
t143 = t147 ^ 2;
t134 = t159 * t144 * t155;
t185 = qJDD(4) + t134;
t163 = t185 * pkin(4) + t199;
t130 = qJD(4) * pkin(4) - pkin(8) * t190;
t150 = t159 ^ 2;
t139 = t150 * t144;
t43 = -pkin(4) * t139 + pkin(8) * t169 - qJD(4) * t130 + t50;
t22 = t154 * t43 - t158 * t163;
t195 = t158 * t43;
t23 = t154 * t163 + t195;
t8 = t154 * t23 - t158 * t22;
t198 = t155 * t8;
t180 = t152 * t90 - t153 * t165;
t56 = -t146 * pkin(3) - t144 * pkin(7) + t180;
t45 = -pkin(4) * t169 - pkin(8) * t139 + t130 * t190 + t56;
t197 = t154 * t45;
t80 = t86 + t145;
t196 = t154 * t80;
t194 = t158 * t45;
t193 = t158 * t80;
t192 = t147 * t154;
t191 = t147 * t158;
t189 = t155 * t185;
t129 = qJDD(4) - t134;
t187 = t159 * t129;
t17 = t152 * t29 - t153 * t56;
t184 = pkin(2) * t17 - pkin(3) * t56 + pkin(7) * t29;
t122 = -t153 * t144 - t152 * t146;
t183 = pkin(2) * t122 - t67;
t9 = t154 * t22 + t158 * t23;
t149 = t155 ^ 2;
t138 = t149 * t144;
t162 = qJD(4) ^ 2;
t131 = -t138 - t162;
t100 = -t155 * t131 - t187;
t118 = 0.2e1 * t181 + t188;
t77 = t152 * t100 - t153 * t118;
t179 = pkin(2) * t77 - pkin(3) * t118 + pkin(7) * t100 + t155 * t56;
t120 = t137 - 0.2e1 * t182;
t132 = -t139 - t162;
t99 = t159 * t132 - t189;
t76 = t153 * t120 + t152 * t99;
t178 = pkin(2) * t76 + pkin(3) * t120 + pkin(7) * t99 - t159 * t56;
t177 = t154 * t119 - t158 * t169;
t164 = (-qJD(5) + t147) * t110 - t177;
t65 = t105 + t72;
t34 = t154 * t164 - t158 * t65;
t35 = t154 * t65 + t158 * t164;
t14 = -t155 * t34 + t159 * t35;
t73 = -t106 - t107;
t11 = t152 * t14 - t153 * t73;
t176 = t155 * (-pkin(8) * t34 - t8) + t159 * (-pkin(4) * t73 + pkin(8) * t35 + t9) - pkin(3) * t73 + pkin(7) * t14 + pkin(2) * t11;
t78 = -t143 - t106;
t51 = t154 * t78 + t202;
t52 = t158 * t78 - t203;
t31 = -t155 * t51 + t159 * t52;
t60 = (qJD(5) + t147) * t110 + t177;
t21 = t152 * t31 - t153 * t60;
t175 = t155 * (-pkin(8) * t51 + t197) + t159 * (-pkin(4) * t60 + pkin(8) * t52 - t194) - pkin(3) * t60 + pkin(7) * t31 + pkin(2) * t21;
t101 = -t107 - t143;
t68 = t158 * t101 - t196;
t69 = -t154 * t101 - t193;
t39 = -t155 * t68 + t159 * t69;
t25 = t152 * t39 - t153 * t201;
t174 = t155 * (-pkin(8) * t68 + t194) + t159 * (-pkin(4) * t201 + pkin(8) * t69 + t197) - pkin(3) * t201 + pkin(7) * t39 + pkin(2) * t25;
t167 = t152 * t144 - t153 * t146;
t173 = -pkin(2) * t167 - t180;
t124 = (t149 + t150) * t146;
t125 = t138 + t139;
t89 = t152 * t124 + t153 * t125;
t172 = pkin(2) * t89 + pkin(3) * t125 + pkin(7) * t124 + t29;
t4 = t159 * t9 - t198;
t2 = t152 * t4 - t153 * t45;
t166 = pkin(2) * t2 + pkin(7) * t4 - pkin(8) * t198 - pkin(3) * t45 + t159 * (-pkin(4) * t45 + pkin(8) * t9);
t103 = -t107 + t143;
t102 = t106 - t143;
t98 = t189 + t159 * (-t138 + t162);
t97 = t155 * (t139 - t162) + t187;
t94 = (t119 + t181) * t155;
t93 = t120 * t159;
t87 = t159 * t118 + t155 * t120;
t82 = t107 - t106;
t71 = -t110 * qJD(5) - t177;
t46 = (t155 * (-t108 * t158 + t110 * t154) + t159 * (-t108 * t154 - t110 * t158)) * t147;
t41 = t155 * (t158 * t102 - t196) + t159 * (t154 * t102 + t193);
t40 = t155 * (-t154 * t103 + t202) + t159 * (t158 * t103 + t203);
t37 = t152 * t67 - t153 * t180;
t36 = pkin(2) * t37;
t33 = t155 * (-t110 * t192 + t158 * t72) + t159 * (t110 * t191 + t154 * t72);
t32 = t155 * (t108 * t191 - t154 * t71) + t159 * (t108 * t192 + t158 * t71);
t13 = t155 * (-t154 * t201 - t158 * t60) + t159 * (-t154 * t60 + t158 * t201);
t1 = [0, 0, 0, 0, 0, qJDD(1), t171, t170, 0, 0, 0, 0, 0, 0, 0, t146, pkin(1) * (-t156 * t144 + t160 * t146) + t91, pkin(1) * (-t160 * t144 - t156 * t146) - t92, 0, pkin(1) * (t156 * t92 + t160 * t91), 0, 0, 0, 0, 0, t146, pkin(1) * (t122 * t156 - t160 * t167) + t173, pkin(1) * (t160 * t122 + t156 * t167) + t183, 0, pkin(1) * (t156 * (t152 * t180 + t153 * t67) + t160 * t37) + t36, t94, t87, t98, t93, t97, 0, pkin(1) * (t156 * (-t152 * t120 + t153 * t99) + t160 * t76) + t178, pkin(1) * (t156 * (t153 * t100 + t152 * t118) + t160 * t77) + t179, pkin(1) * (t156 * (t153 * t124 - t152 * t125) + t160 * t89) + t172, pkin(1) * (t156 * (t152 * t56 + t153 * t29) + t160 * t17) + t184, t33, t13, t40, t32, t41, t46, pkin(1) * (t156 * (t152 * t60 + t153 * t31) + t160 * t21) + t175, pkin(1) * (t156 * (t152 * t201 + t153 * t39) + t160 * t25) + t174, pkin(1) * (t156 * (t153 * t14 + t152 * t73) + t160 * t11) + t176, pkin(1) * (t156 * (t152 * t45 + t153 * t4) + t160 * t2) + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t91, -t92, 0, 0, 0, 0, 0, 0, 0, t146, t173, t183, 0, t36, t94, t87, t98, t93, t97, 0, t178, t179, t172, t184, t33, t13, t40, t32, t41, t46, t175, t174, t176, t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, 0, 0, 0, 0, 0, 0, t155 * t132 + t159 * t185, -t155 * t129 + t159 * t131, 0, t155 * t50 - t159 * t49, 0, 0, 0, 0, 0, 0, t155 * t52 + t159 * t51, t155 * t69 + t159 * t68, t155 * t35 + t159 * t34, t155 * t9 + t159 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, t138 - t139, t188, t134, t137, qJDD(4), -t49, -t50, 0, 0, t86, t82, t65, -t86, t164, t145, pkin(4) * t51 - t22, -t195 - t154 * t199 + (-t154 * t185 + t68) * pkin(4), pkin(4) * t34, pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t82, t65, -t86, t164, t145, -t22, -t23, 0, 0;];
tauJ_reg = t1;
