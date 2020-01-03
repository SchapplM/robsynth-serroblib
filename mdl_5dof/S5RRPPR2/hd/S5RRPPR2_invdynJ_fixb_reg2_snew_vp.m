% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:36
% EndTime: 2020-01-03 11:57:40
% DurationCPUTime: 1.09s
% Computational Cost: add. (4826->165), mult. (6694->245), div. (0->0), fcn. (3931->10), ass. (0->123)
t128 = sin(pkin(9));
t132 = sin(qJ(5));
t125 = qJD(1) + qJD(2);
t162 = qJD(5) * t125;
t106 = t132 * t128 * t162;
t130 = cos(pkin(9));
t167 = t125 * t130;
t110 = -qJD(5) + t167;
t155 = t132 * t110 * t125;
t122 = qJDD(1) + qJDD(2);
t135 = cos(qJ(5));
t168 = t122 * t135;
t140 = (t155 + t168) * t128 - t106;
t179 = t140 * t128;
t173 = t130 * pkin(4);
t174 = t128 * pkin(7);
t148 = -t173 - t174;
t175 = 2 * qJD(4);
t121 = t125 ^ 2;
t129 = sin(pkin(8));
t131 = cos(pkin(8));
t134 = sin(qJ(1));
t137 = cos(qJ(1));
t147 = -g(2) * t137 - t134 * g(3);
t103 = qJDD(1) * pkin(1) + t147;
t146 = t134 * g(2) - g(3) * t137;
t104 = -qJD(1) ^ 2 * pkin(1) - t146;
t133 = sin(qJ(2));
t136 = cos(qJ(2));
t78 = t136 * t103 - t133 * t104;
t72 = pkin(2) * t122 + t78;
t79 = t133 * t103 + t136 * t104;
t73 = -pkin(2) * t121 + t79;
t54 = t129 * t72 + t131 * t73;
t49 = -pkin(3) * t121 + qJ(4) * t122 + t54;
t178 = t125 * t175 + t49;
t123 = t128 ^ 2;
t124 = t130 ^ 2;
t100 = (t123 + t124) * t121;
t163 = -g(1) + qJDD(3);
t118 = t130 * t163;
t35 = t178 * t128 - t118;
t36 = t128 * t163 + t178 * t130;
t16 = t128 * t35 + t130 * t36;
t169 = t122 * t132;
t177 = (t135 * t162 + t169) * t128;
t107 = t110 ^ 2;
t170 = t121 * t123;
t156 = t132 * t170;
t102 = t135 * t156;
t164 = t130 * t122;
t108 = -qJDD(5) + t164;
t81 = -t102 + t108;
t172 = t132 * t81;
t82 = -t102 - t108;
t171 = t135 * t82;
t166 = t125 * t135;
t165 = t128 * t122;
t94 = t148 * t125;
t27 = -t118 + (t49 + (t175 + t94) * t125) * t128;
t28 = t94 * t167 + t36;
t154 = t129 * t73 - t131 * t72;
t48 = -t122 * pkin(3) - t121 * qJ(4) + qJDD(4) + t154;
t37 = t148 * t122 + t48;
t13 = t132 * t28 - t135 * t37;
t14 = t132 * t37 + t135 * t28;
t8 = t13 * t132 + t135 * t14;
t4 = t128 * t27 + t130 * t8;
t7 = -t13 * t135 + t132 * t14;
t2 = t129 * t4 - t131 * t7;
t161 = pkin(2) * t2 - pkin(3) * t7 + qJ(4) * t4;
t10 = t129 * t16 - t131 * t48;
t160 = pkin(2) * t10 - pkin(3) * t48 + qJ(4) * t16;
t98 = -t121 * t131 - t122 * t129;
t159 = pkin(2) * t98 - t54;
t157 = t135 ^ 2 * t170;
t90 = t130 * t100;
t77 = -t129 * t90 + t131 * t164;
t153 = pkin(2) * t77 + pkin(3) * t164 - qJ(4) * t90 - t130 * t48;
t145 = t121 * t129 - t122 * t131;
t152 = -pkin(2) * t145 - t154;
t80 = -t157 - t107;
t59 = -t132 * t80 + t135 * t81;
t39 = t130 * t59 + t179;
t58 = t135 * t80 + t172;
t22 = t129 * t39 - t131 * t58;
t151 = t128 * (-pkin(7) * t58 + t135 * t27) + t130 * (-pkin(4) * t58 + t14) - pkin(3) * t58 + qJ(4) * t39 + pkin(2) * t22;
t109 = t132 ^ 2 * t170;
t83 = -t109 - t107;
t62 = -t132 * t82 + t135 * t83;
t88 = t110 * t128 * t166;
t64 = t88 - t177;
t43 = -t128 * t64 + t130 * t62;
t61 = t132 * t83 + t171;
t26 = t129 * t43 - t131 * t61;
t150 = t128 * (-pkin(7) * t61 + t132 * t27) + t130 * (-pkin(4) * t61 + t13) - pkin(3) * t61 + qJ(4) * t43 + pkin(2) * t26;
t115 = t123 * t122;
t116 = t124 * t122;
t96 = t116 + t115;
t71 = t100 * t131 + t129 * t96;
t149 = pkin(2) * t71 + pkin(3) * t100 + qJ(4) * t96 + t16;
t89 = t128 * t100;
t76 = t129 * t89 - t131 * t165;
t144 = pkin(2) * t76 - pkin(3) * t165 + qJ(4) * t89 + t128 * t48;
t63 = t88 + t177;
t65 = -t106 + (-t155 + t168) * t128;
t52 = t132 * t65 - t135 * t63;
t84 = t109 + t157;
t31 = -t128 * t84 + t130 * t52;
t51 = -t132 * t63 - t135 * t65;
t18 = t129 * t31 - t131 * t51;
t141 = pkin(2) * t18 + qJ(4) * t31 - t128 * t7 + (-pkin(3) + t148) * t51;
t105 = 0.2e1 * t128 * t164;
t101 = t130 * t108;
t85 = -t109 + t157;
t56 = (-t130 * t156 + t179) * t135;
t55 = t130 * t102 + (t169 + (qJD(5) - t110) * t166) * t132 * t123;
t42 = t128 * (t135 * (t109 - t107) + t172) + t130 * t63;
t41 = t128 * (t171 - t132 * (t107 - t157)) - t130 * t65;
t30 = t128 * (-t132 * t140 + t135 * t64) - t130 * t85;
t24 = t129 * t54 - t131 * t154;
t23 = pkin(2) * t24;
t1 = [0, 0, 0, 0, 0, qJDD(1), t147, t146, 0, 0, 0, 0, 0, 0, 0, t122, pkin(1) * (-t121 * t133 + t122 * t136) + t78, pkin(1) * (-t121 * t136 - t122 * t133) - t79, 0, pkin(1) * (t133 * t79 + t136 * t78), 0, 0, 0, 0, 0, t122, pkin(1) * (t133 * t98 - t136 * t145) + t152, pkin(1) * (t133 * t145 + t136 * t98) + t159, 0, pkin(1) * (t133 * (t129 * t154 + t131 * t54) + t136 * t24) + t23, t115, t105, 0, t116, 0, 0, pkin(1) * (t133 * (-t129 * t164 - t131 * t90) + t136 * t77) + t153, pkin(1) * (t133 * (t129 * t165 + t131 * t89) + t136 * t76) + t144, pkin(1) * (t133 * (-t100 * t129 + t131 * t96) + t136 * t71) + t149, pkin(1) * (t133 * (t129 * t48 + t131 * t16) + t136 * t10) + t160, t56, t30, t41, t55, t42, t101, pkin(1) * (t133 * (t129 * t61 + t131 * t43) + t136 * t26) + t150, pkin(1) * (t133 * (t129 * t58 + t131 * t39) + t136 * t22) + t151, pkin(1) * (t133 * (t129 * t51 + t131 * t31) + t136 * t18) + t141, pkin(1) * (t133 * (t129 * t7 + t131 * t4) + t136 * t2) - t7 * t174 - t7 * t173 + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t78, -t79, 0, 0, 0, 0, 0, 0, 0, t122, t152, t159, 0, t23, t115, t105, 0, t116, 0, 0, t153, t144, t149, t160, t56, t30, t41, t55, t42, t101, t150, t151, t141, t148 * t7 + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128 * t36 - t130 * t35, 0, 0, 0, 0, 0, 0, t128 * t62 + t130 * t64, t128 * t59 - t130 * t140, t128 * t52 + t130 * t84, t128 * t8 - t130 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t165, -t100, t48, 0, 0, 0, 0, 0, 0, t61, t58, t51, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t85, t65, -t102, -t63, -t108, -t13, -t14, 0, 0;];
tauJ_reg = t1;
