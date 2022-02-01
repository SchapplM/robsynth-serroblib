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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:05:52
% EndTime: 2022-01-20 10:05:56
% DurationCPUTime: 1.14s
% Computational Cost: add. (4826->165), mult. (6694->245), div. (0->0), fcn. (3931->10), ass. (0->123)
t129 = sin(pkin(9));
t133 = sin(qJ(5));
t126 = qJD(1) + qJD(2);
t163 = qJD(5) * t126;
t106 = t133 * t129 * t163;
t131 = cos(pkin(9));
t168 = t126 * t131;
t110 = -qJD(5) + t168;
t158 = t110 * t126 * t133;
t123 = qJDD(1) + qJDD(2);
t136 = cos(qJ(5));
t169 = t123 * t136;
t141 = (t158 + t169) * t129 - t106;
t180 = t141 * t129;
t174 = pkin(7) * t129;
t175 = pkin(4) * t131;
t148 = -t174 - t175;
t176 = 2 * qJD(4);
t122 = t126 ^ 2;
t130 = sin(pkin(8));
t132 = cos(pkin(8));
t135 = sin(qJ(1));
t138 = cos(qJ(1));
t155 = t135 * g(1) - g(2) * t138;
t103 = qJDD(1) * pkin(1) + t155;
t147 = g(1) * t138 + t135 * g(2);
t104 = -qJD(1) ^ 2 * pkin(1) - t147;
t134 = sin(qJ(2));
t137 = cos(qJ(2));
t78 = t137 * t103 - t134 * t104;
t72 = pkin(2) * t123 + t78;
t79 = t134 * t103 + t137 * t104;
t73 = -pkin(2) * t122 + t79;
t54 = t130 * t72 + t132 * t73;
t49 = -pkin(3) * t122 + qJ(4) * t123 + t54;
t179 = t126 * t176 + t49;
t124 = t129 ^ 2;
t125 = t131 ^ 2;
t100 = (t124 + t125) * t122;
t164 = -g(3) + qJDD(3);
t118 = t131 * t164;
t35 = t179 * t129 - t118;
t36 = t129 * t164 + t179 * t131;
t16 = t129 * t35 + t131 * t36;
t170 = t123 * t133;
t178 = (t136 * t163 + t170) * t129;
t107 = t110 ^ 2;
t171 = t122 * t124;
t156 = t133 * t171;
t102 = t136 * t156;
t165 = t131 * t123;
t108 = -qJDD(5) + t165;
t81 = -t102 + t108;
t173 = t133 * t81;
t82 = -t102 - t108;
t172 = t136 * t82;
t167 = t126 * t136;
t166 = t129 * t123;
t94 = t148 * t126;
t27 = -t118 + (t49 + (t176 + t94) * t126) * t129;
t28 = t94 * t168 + t36;
t154 = t130 * t73 - t132 * t72;
t48 = -t123 * pkin(3) - t122 * qJ(4) + qJDD(4) + t154;
t37 = t148 * t123 + t48;
t13 = t133 * t28 - t136 * t37;
t14 = t133 * t37 + t136 * t28;
t8 = t13 * t133 + t136 * t14;
t4 = t129 * t27 + t131 * t8;
t7 = -t13 * t136 + t133 * t14;
t2 = t130 * t4 - t132 * t7;
t162 = pkin(2) * t2 - pkin(3) * t7 + qJ(4) * t4;
t10 = t130 * t16 - t132 * t48;
t161 = pkin(2) * t10 - pkin(3) * t48 + qJ(4) * t16;
t98 = -t122 * t132 - t123 * t130;
t160 = pkin(2) * t98 - t54;
t157 = t136 ^ 2 * t171;
t90 = t131 * t100;
t77 = -t130 * t90 + t132 * t165;
t153 = pkin(2) * t77 + pkin(3) * t165 - qJ(4) * t90 - t131 * t48;
t146 = t122 * t130 - t123 * t132;
t152 = -pkin(2) * t146 - t154;
t80 = -t157 - t107;
t59 = -t133 * t80 + t136 * t81;
t39 = t131 * t59 + t180;
t58 = t136 * t80 + t173;
t22 = t130 * t39 - t132 * t58;
t151 = t129 * (-pkin(7) * t58 + t136 * t27) + t131 * (-pkin(4) * t58 + t14) - pkin(3) * t58 + qJ(4) * t39 + pkin(2) * t22;
t109 = t133 ^ 2 * t171;
t83 = -t109 - t107;
t62 = -t133 * t82 + t136 * t83;
t88 = t110 * t129 * t167;
t64 = t88 - t178;
t43 = -t129 * t64 + t131 * t62;
t61 = t133 * t83 + t172;
t26 = t130 * t43 - t132 * t61;
t150 = t129 * (-pkin(7) * t61 + t133 * t27) + t131 * (-pkin(4) * t61 + t13) - pkin(3) * t61 + qJ(4) * t43 + pkin(2) * t26;
t115 = t124 * t123;
t116 = t125 * t123;
t96 = t116 + t115;
t71 = t100 * t132 + t130 * t96;
t149 = pkin(2) * t71 + pkin(3) * t100 + qJ(4) * t96 + t16;
t89 = t129 * t100;
t76 = t130 * t89 - t132 * t166;
t145 = pkin(2) * t76 - pkin(3) * t166 + qJ(4) * t89 + t129 * t48;
t63 = t88 + t178;
t65 = -t106 + (-t158 + t169) * t129;
t52 = t133 * t65 - t136 * t63;
t84 = t109 + t157;
t31 = -t129 * t84 + t131 * t52;
t51 = -t133 * t63 - t136 * t65;
t18 = t130 * t31 - t132 * t51;
t142 = pkin(2) * t18 + qJ(4) * t31 - t129 * t7 + (-pkin(3) + t148) * t51;
t105 = 0.2e1 * t129 * t165;
t101 = t131 * t108;
t85 = -t109 + t157;
t56 = (-t131 * t156 + t180) * t136;
t55 = t131 * t102 + (t170 + (qJD(5) - t110) * t167) * t133 * t124;
t42 = t129 * (t136 * (t109 - t107) + t173) + t131 * t63;
t41 = t129 * (t172 - t133 * (t107 - t157)) - t131 * t65;
t30 = t129 * (-t133 * t141 + t136 * t64) - t131 * t85;
t24 = t130 * t54 - t132 * t154;
t23 = pkin(2) * t24;
t1 = [0, 0, 0, 0, 0, qJDD(1), t155, t147, 0, 0, 0, 0, 0, 0, 0, t123, pkin(1) * (-t122 * t134 + t123 * t137) + t78, pkin(1) * (-t122 * t137 - t123 * t134) - t79, 0, pkin(1) * (t134 * t79 + t137 * t78), 0, 0, 0, 0, 0, t123, pkin(1) * (t134 * t98 - t137 * t146) + t152, pkin(1) * (t134 * t146 + t137 * t98) + t160, 0, pkin(1) * (t134 * (t130 * t154 + t132 * t54) + t137 * t24) + t23, t115, t105, 0, t116, 0, 0, pkin(1) * (t134 * (-t130 * t165 - t132 * t90) + t137 * t77) + t153, pkin(1) * (t134 * (t130 * t166 + t132 * t89) + t137 * t76) + t145, pkin(1) * (t134 * (-t100 * t130 + t132 * t96) + t137 * t71) + t149, pkin(1) * (t134 * (t130 * t48 + t132 * t16) + t137 * t10) + t161, t56, t30, t41, t55, t42, t101, pkin(1) * (t134 * (t130 * t61 + t132 * t43) + t137 * t26) + t150, pkin(1) * (t134 * (t130 * t58 + t132 * t39) + t137 * t22) + t151, pkin(1) * (t134 * (t130 * t51 + t132 * t31) + t137 * t18) + t142, pkin(1) * (t134 * (t130 * t7 + t132 * t4) + t137 * t2) - t7 * t174 - t7 * t175 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t78, -t79, 0, 0, 0, 0, 0, 0, 0, t123, t152, t160, 0, t23, t115, t105, 0, t116, 0, 0, t153, t145, t149, t161, t56, t30, t41, t55, t42, t101, t150, t151, t142, t148 * t7 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129 * t36 - t131 * t35, 0, 0, 0, 0, 0, 0, t129 * t62 + t131 * t64, t129 * t59 - t131 * t141, t129 * t52 + t131 * t84, t129 * t8 - t131 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, t166, -t100, t48, 0, 0, 0, 0, 0, 0, t61, t58, t51, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t85, t65, -t102, -t63, -t108, -t13, -t14, 0, 0;];
tauJ_reg = t1;
