% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:40
% EndTime: 2019-12-31 20:15:44
% DurationCPUTime: 1.18s
% Computational Cost: add. (5019->179), mult. (6254->235), div. (0->0), fcn. (3771->8), ass. (0->136)
t129 = sin(qJ(5));
t123 = qJDD(4) + qJDD(5);
t126 = qJD(1) + qJD(2);
t130 = sin(qJ(4));
t133 = cos(qJ(5));
t134 = cos(qJ(4));
t86 = (-t129 * t134 - t130 * t133) * t126;
t169 = t126 * t134;
t88 = -t129 * t130 * t126 + t133 * t169;
t184 = t88 * t86;
t192 = t123 + t184;
t196 = t129 * t192;
t195 = t133 * t192;
t194 = pkin(2) + pkin(7);
t125 = qJD(4) + qJD(5);
t182 = t125 * t86;
t164 = qJD(4) * t126;
t161 = t134 * t164;
t124 = qJDD(1) + qJDD(2);
t167 = t130 * t124;
t98 = -t161 - t167;
t112 = t134 * t124;
t160 = t130 * t164;
t99 = t112 - t160;
t49 = t86 * qJD(5) + t129 * t98 + t133 * t99;
t193 = t49 + t182;
t122 = t126 ^ 2;
t131 = sin(qJ(2));
t135 = cos(qJ(2));
t191 = pkin(1) * (t122 * t131 - t124 * t135);
t190 = pkin(1) * (t122 * t135 + t124 * t131);
t186 = t124 * pkin(2);
t132 = sin(qJ(1));
t136 = cos(qJ(1));
t159 = t132 * g(1) - g(2) * t136;
t105 = qJDD(1) * pkin(1) + t159;
t155 = g(1) * t136 + g(2) * t132;
t106 = -qJD(1) ^ 2 * pkin(1) - t155;
t69 = t135 * t105 - t131 * t106;
t156 = qJDD(3) - t186 - t69;
t62 = -t122 * qJ(3) + t156;
t144 = -pkin(7) * t124 + t62;
t139 = t134 * t144;
t52 = t130 * g(3) + t139;
t53 = g(3) * t134 - t130 * t144;
t29 = -t130 * t53 + t134 * t52;
t173 = t122 * t134;
t187 = t99 * pkin(8);
t189 = qJDD(4) * pkin(4) + t139 + (-pkin(4) * t173 - pkin(8) * t164 + g(3)) * t130 - t187;
t84 = t86 ^ 2;
t85 = t88 ^ 2;
t121 = t125 ^ 2;
t188 = 2 * qJD(3);
t109 = qJD(4) * pkin(4) - pkin(8) * t169;
t127 = t130 ^ 2;
t175 = t122 * t127;
t34 = -pkin(4) * t175 + pkin(8) * t98 - qJD(4) * t109 - t53;
t14 = t129 * t34 - t133 * t189;
t177 = t133 * t34;
t15 = t189 * t129 + t177;
t6 = t129 * t15 - t133 * t14;
t185 = t134 * t6;
t172 = t124 * qJ(3);
t70 = t131 * t105 + t135 * t106;
t143 = -t122 * pkin(2) + t172 + t70;
t163 = t126 * t188;
t59 = t143 + t163;
t183 = -pkin(2) * t62 + qJ(3) * t59;
t118 = t122 * pkin(7);
t33 = -pkin(8) * t175 - t118 - t98 * pkin(4) + (t109 * t134 + t188) * t126 + t143;
t181 = t129 * t33;
t64 = -t184 + t123;
t180 = t129 * t64;
t178 = t133 * t33;
t176 = t133 * t64;
t128 = t134 ^ 2;
t174 = t122 * t128;
t171 = t125 * t129;
t170 = t125 * t133;
t162 = t130 * t173;
t168 = t130 * (qJDD(4) + t162);
t108 = qJDD(4) - t162;
t166 = t134 * t108;
t165 = t127 + t128;
t7 = t129 * t14 + t133 * t15;
t158 = t129 * t99 - t133 * t98;
t56 = -t118 + t59;
t157 = qJ(3) * t56 - t194 * t29;
t100 = t112 - 0.2e1 * t160;
t137 = qJD(4) ^ 2;
t76 = -t168 + t134 * (-t137 - t174);
t154 = qJ(3) * t100 + t134 * t56 - t194 * t76;
t151 = -t186 + t156;
t75 = t130 * (-t137 - t175) + t166;
t97 = 0.2e1 * t161 + t167;
t150 = qJ(3) * t97 + t130 * t56 - t194 * t75;
t145 = (-qJD(5) + t125) * t88 - t158;
t44 = -t182 + t49;
t22 = t129 * t145 - t133 * t44;
t23 = t129 * t44 + t133 * t145;
t51 = -t84 - t85;
t9 = t130 * t23 + t134 * t22;
t149 = -t130 * (-pkin(4) * t51 + pkin(8) * t23 + t7) + t134 * (-pkin(8) * t22 - t6) + qJ(3) * t51 - t194 * t9;
t61 = -t121 - t84;
t35 = t129 * t61 + t195;
t36 = t133 * t61 - t196;
t18 = t130 * t36 + t134 * t35;
t39 = (qJD(5) + t125) * t88 + t158;
t148 = -t130 * (-pkin(4) * t39 + pkin(8) * t36 - t178) + t134 * (-pkin(8) * t35 + t181) + qJ(3) * t39 - t194 * t18;
t79 = -t85 - t121;
t45 = t133 * t79 - t180;
t46 = -t129 * t79 - t176;
t25 = t130 * t46 + t134 * t45;
t147 = -t130 * (-pkin(4) * t193 + pkin(8) * t46 + t181) + t134 * (-pkin(8) * t45 + t178) + qJ(3) * t193 - t194 * t25;
t103 = t165 * t124;
t104 = t165 * t122;
t146 = -qJ(3) * t104 + t194 * t103 - t29;
t142 = 0.2e1 * t172 + t70 + t163;
t2 = t130 * t7 + t185;
t141 = -pkin(8) * t185 - t130 * (-pkin(4) * t33 + pkin(8) * t7) + qJ(3) * t33 - t194 * t2;
t81 = -t85 + t121;
t80 = t84 - t121;
t78 = t166 - t130 * (t137 - t174);
t77 = t134 * (-t137 + t175) - t168;
t72 = (t99 - t160) * t134;
t71 = (-t98 + t161) * t130;
t68 = -t100 * t130 - t134 * t97;
t66 = t85 - t84;
t48 = -qJD(5) * t88 - t158;
t30 = (t134 * (t129 * t88 + t133 * t86) - t130 * (t129 * t86 - t133 * t88)) * t125;
t27 = t134 * (t133 * t80 - t180) - t130 * (t129 * t80 + t176);
t26 = t134 * (-t129 * t81 + t195) - t130 * (t133 * t81 + t196);
t20 = t134 * (t133 * t49 - t88 * t171) - t130 * (t129 * t49 + t88 * t170);
t19 = t134 * (-t129 * t48 - t86 * t170) - t130 * (t133 * t48 - t86 * t171);
t10 = t134 * (-t129 * t193 - t133 * t39) - t130 * (-t129 * t39 + t133 * t193);
t1 = [0, 0, 0, 0, 0, qJDD(1), t159, t155, 0, 0, 0, 0, 0, 0, 0, t124, t69 - t191, -t70 - t190, 0, pkin(1) * (t131 * t70 + t135 * t69), t124, 0, 0, 0, 0, 0, 0, t151 + t191, t142 + t190, pkin(1) * (t131 * t59 - t135 * t62) + t183, t72, t68, t78, t71, t77, 0, pkin(1) * (t131 * t97 - t135 * t75) + t150, pkin(1) * (t100 * t131 - t135 * t76) + t154, pkin(1) * (t103 * t135 - t104 * t131) + t146, pkin(1) * (t131 * t56 - t135 * t29) + t157, t20, t10, t26, t19, t27, t30, pkin(1) * (t131 * t39 - t135 * t18) + t148, pkin(1) * (t131 * t193 - t135 * t25) + t147, pkin(1) * (t131 * t51 - t135 * t9) + t149, pkin(1) * (t131 * t33 - t135 * t2) + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t69, -t70, 0, 0, t124, 0, 0, 0, 0, 0, 0, t151, t142, t183, t72, t68, t78, t71, t77, 0, t150, t154, t146, t157, t20, t10, t26, t19, t27, t30, t148, t147, t149, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, -t122, t62, 0, 0, 0, 0, 0, 0, t75, t76, -t103, t29, 0, 0, 0, 0, 0, 0, t18, t25, t9, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, (-t127 + t128) * t122, t112, -t162, -t167, qJDD(4), t52, t53, 0, 0, -t184, t66, t44, t184, t145, t123, pkin(4) * t35 - t14, -t177 - t129 * (-pkin(8) * t160 - t187 + t52) + (-t129 * t108 + t45) * pkin(4), pkin(4) * t22, pkin(4) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, t66, t44, t184, t145, t123, -t14, -t15, 0, 0;];
tauJ_reg = t1;
