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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:20:22
% EndTime: 2019-12-05 18:20:28
% DurationCPUTime: 1.06s
% Computational Cost: add. (4826->165), mult. (6694->245), div. (0->0), fcn. (3931->10), ass. (0->123)
t130 = sin(pkin(9));
t134 = sin(qJ(5));
t127 = qJD(1) + qJD(2);
t163 = qJD(5) * t127;
t106 = t134 * t130 * t163;
t132 = cos(pkin(9));
t169 = t127 * t132;
t110 = -qJD(5) + t169;
t158 = t110 * t127 * t134;
t124 = qJDD(1) + qJDD(2);
t137 = cos(qJ(5));
t170 = t124 * t137;
t142 = (t158 + t170) * t130 - t106;
t181 = t142 * t130;
t175 = t132 * pkin(4);
t176 = t130 * pkin(7);
t149 = -t175 - t176;
t177 = 2 * qJD(4);
t123 = t127 ^ 2;
t131 = sin(pkin(8));
t133 = cos(pkin(8));
t136 = sin(qJ(1));
t139 = cos(qJ(1));
t164 = t139 * g(2) + t136 * g(3);
t103 = qJDD(1) * pkin(1) + t164;
t148 = -t136 * g(2) + g(3) * t139;
t104 = -qJD(1) ^ 2 * pkin(1) - t148;
t135 = sin(qJ(2));
t138 = cos(qJ(2));
t78 = t138 * t103 - t135 * t104;
t72 = pkin(2) * t124 + t78;
t79 = t135 * t103 + t138 * t104;
t73 = -pkin(2) * t123 + t79;
t54 = t131 * t72 + t133 * t73;
t49 = -pkin(3) * t123 + qJ(4) * t124 + t54;
t180 = t127 * t177 + t49;
t125 = t130 ^ 2;
t126 = t132 ^ 2;
t100 = (t125 + t126) * t123;
t165 = -g(1) + qJDD(3);
t118 = t132 * t165;
t35 = t180 * t130 - t118;
t36 = t130 * t165 + t180 * t132;
t16 = t130 * t35 + t132 * t36;
t171 = t124 * t134;
t179 = (t137 * t163 + t171) * t130;
t107 = t110 ^ 2;
t172 = t123 * t125;
t156 = t134 * t172;
t102 = t137 * t156;
t166 = t132 * t124;
t108 = -qJDD(5) + t166;
t81 = -t102 + t108;
t174 = t134 * t81;
t82 = -t102 - t108;
t173 = t137 * t82;
t168 = t127 * t137;
t167 = t130 * t124;
t94 = t149 * t127;
t27 = -t118 + (t49 + (t177 + t94) * t127) * t130;
t28 = t169 * t94 + t36;
t155 = t131 * t73 - t133 * t72;
t48 = -t124 * pkin(3) - t123 * qJ(4) + qJDD(4) + t155;
t37 = t124 * t149 + t48;
t13 = t134 * t28 - t137 * t37;
t14 = t134 * t37 + t137 * t28;
t8 = t13 * t134 + t137 * t14;
t4 = t130 * t27 + t132 * t8;
t7 = -t13 * t137 + t134 * t14;
t2 = t131 * t4 - t133 * t7;
t162 = pkin(2) * t2 - pkin(3) * t7 + qJ(4) * t4;
t10 = t131 * t16 - t133 * t48;
t161 = pkin(2) * t10 - pkin(3) * t48 + qJ(4) * t16;
t98 = -t123 * t133 - t124 * t131;
t160 = pkin(2) * t98 - t54;
t157 = t137 ^ 2 * t172;
t90 = t132 * t100;
t77 = -t131 * t90 + t133 * t166;
t154 = pkin(2) * t77 + pkin(3) * t166 - qJ(4) * t90 - t132 * t48;
t147 = t123 * t131 - t124 * t133;
t153 = -pkin(2) * t147 - t155;
t80 = -t157 - t107;
t59 = -t134 * t80 + t137 * t81;
t39 = t132 * t59 + t181;
t58 = t137 * t80 + t174;
t22 = t131 * t39 - t133 * t58;
t152 = t130 * (-pkin(7) * t58 + t137 * t27) + t132 * (-pkin(4) * t58 + t14) - pkin(3) * t58 + qJ(4) * t39 + pkin(2) * t22;
t109 = t134 ^ 2 * t172;
t83 = -t109 - t107;
t62 = -t134 * t82 + t137 * t83;
t88 = t110 * t130 * t168;
t64 = t88 - t179;
t43 = -t130 * t64 + t132 * t62;
t61 = t134 * t83 + t173;
t26 = t131 * t43 - t133 * t61;
t151 = t130 * (-pkin(7) * t61 + t134 * t27) + t132 * (-pkin(4) * t61 + t13) - pkin(3) * t61 + qJ(4) * t43 + pkin(2) * t26;
t115 = t125 * t124;
t116 = t126 * t124;
t96 = t116 + t115;
t71 = t100 * t133 + t131 * t96;
t150 = pkin(2) * t71 + pkin(3) * t100 + qJ(4) * t96 + t16;
t89 = t130 * t100;
t76 = t131 * t89 - t133 * t167;
t146 = pkin(2) * t76 - pkin(3) * t167 + qJ(4) * t89 + t130 * t48;
t63 = t88 + t179;
t65 = -t106 + (-t158 + t170) * t130;
t52 = t134 * t65 - t137 * t63;
t84 = t109 + t157;
t31 = -t130 * t84 + t132 * t52;
t51 = -t134 * t63 - t137 * t65;
t18 = t131 * t31 - t133 * t51;
t143 = pkin(2) * t18 + qJ(4) * t31 - t130 * t7 + (-pkin(3) + t149) * t51;
t105 = 0.2e1 * t130 * t166;
t101 = t132 * t108;
t85 = -t109 + t157;
t56 = (-t132 * t156 + t181) * t137;
t55 = t132 * t102 + (t171 + (qJD(5) - t110) * t168) * t134 * t125;
t42 = t130 * (t137 * (t109 - t107) + t174) + t132 * t63;
t41 = t130 * (t173 - t134 * (t107 - t157)) - t132 * t65;
t30 = t130 * (-t134 * t142 + t137 * t64) - t132 * t85;
t24 = t131 * t54 - t133 * t155;
t23 = pkin(2) * t24;
t1 = [0, 0, 0, 0, 0, qJDD(1), t164, t148, 0, 0, 0, 0, 0, 0, 0, t124, pkin(1) * (-t123 * t135 + t124 * t138) + t78, pkin(1) * (-t123 * t138 - t124 * t135) - t79, 0, pkin(1) * (t135 * t79 + t138 * t78), 0, 0, 0, 0, 0, t124, pkin(1) * (t135 * t98 - t138 * t147) + t153, pkin(1) * (t135 * t147 + t138 * t98) + t160, 0, pkin(1) * (t135 * (t131 * t155 + t133 * t54) + t138 * t24) + t23, t115, t105, 0, t116, 0, 0, pkin(1) * (t135 * (-t131 * t166 - t133 * t90) + t138 * t77) + t154, pkin(1) * (t135 * (t131 * t167 + t133 * t89) + t138 * t76) + t146, pkin(1) * (t135 * (-t100 * t131 + t133 * t96) + t138 * t71) + t150, pkin(1) * (t135 * (t131 * t48 + t133 * t16) + t138 * t10) + t161, t56, t30, t41, t55, t42, t101, pkin(1) * (t135 * (t131 * t61 + t133 * t43) + t138 * t26) + t151, pkin(1) * (t135 * (t131 * t58 + t133 * t39) + t138 * t22) + t152, pkin(1) * (t135 * (t131 * t51 + t133 * t31) + t138 * t18) + t143, pkin(1) * (t135 * (t131 * t7 + t133 * t4) + t138 * t2) - t7 * t176 - t7 * t175 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t78, -t79, 0, 0, 0, 0, 0, 0, 0, t124, t153, t160, 0, t23, t115, t105, 0, t116, 0, 0, t154, t146, t150, t161, t56, t30, t41, t55, t42, t101, t151, t152, t143, t149 * t7 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130 * t36 - t132 * t35, 0, 0, 0, 0, 0, 0, t130 * t62 + t132 * t64, t130 * t59 - t132 * t142, t130 * t52 + t132 * t84, t130 * t8 - t132 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, t167, -t100, t48, 0, 0, 0, 0, 0, 0, t61, t58, t51, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t85, t65, -t102, -t63, -t108, -t13, -t14, 0, 0;];
tauJ_reg = t1;
