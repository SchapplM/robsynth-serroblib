% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:26
% EndTime: 2019-12-31 21:49:29
% DurationCPUTime: 0.96s
% Computational Cost: add. (3714->150), mult. (4298->178), div. (0->0), fcn. (2296->8), ass. (0->112)
t107 = sin(qJ(3));
t108 = sin(qJ(2));
t111 = cos(qJ(3));
t112 = cos(qJ(2));
t106 = sin(qJ(4));
t110 = cos(qJ(4));
t100 = qJD(1) + qJD(2);
t95 = qJD(3) + t100;
t93 = t95 ^ 2;
t87 = t110 * t93 * t106;
t80 = qJDD(4) - t87;
t149 = t110 * t80;
t114 = qJD(4) ^ 2;
t103 = t106 ^ 2;
t158 = t103 * t93;
t83 = t114 + t158;
t54 = -t106 * t83 + t149;
t146 = qJD(4) * t95;
t99 = qJDD(1) + qJDD(2);
t94 = qJDD(3) + t99;
t153 = t106 * t94;
t68 = 0.2e1 * t110 * t146 + t153;
t37 = t107 * t54 + t111 * t68;
t174 = pkin(1) * (t108 * (-t107 * t68 + t111 * t54) + t112 * t37);
t173 = pkin(2) * t37;
t170 = pkin(3) * t68 + pkin(8) * t54;
t168 = qJ(5) * t68;
t162 = t110 * g(3);
t109 = sin(qJ(1));
t113 = cos(qJ(1));
t138 = t109 * g(1) - t113 * g(2);
t81 = qJDD(1) * pkin(1) + t138;
t128 = t113 * g(1) + t109 * g(2);
t82 = -qJD(1) ^ 2 * pkin(1) - t128;
t50 = -t108 * t82 + t112 * t81;
t125 = t99 * pkin(2) + t50;
t51 = t108 * t81 + t112 * t82;
t98 = t100 ^ 2;
t44 = -t98 * pkin(2) + t51;
t32 = t107 * t125 + t111 * t44;
t29 = -t93 * pkin(3) + t94 * pkin(8) + t32;
t22 = t106 * t29 + t162;
t23 = -t106 * g(3) + t110 * t29;
t8 = t106 * t22 + t110 * t23;
t104 = t110 ^ 2;
t157 = t104 * t93;
t55 = t149 + t106 * (-t114 + t157);
t145 = t106 * qJ(5);
t163 = t110 * pkin(4);
t127 = -t145 - t163;
t161 = t127 * t93;
t18 = -qJDD(4) * pkin(4) - t114 * qJ(5) + (t29 + t161) * t106 + qJDD(5) + t162;
t167 = 2 * qJD(5);
t136 = t107 * t44 - t111 * t125;
t28 = -t94 * pkin(3) - t93 * pkin(8) + t136;
t165 = -pkin(3) * t28 + pkin(8) * t8;
t140 = t106 * t146;
t148 = t110 * t94;
t118 = t28 - (-t140 + t148) * pkin(4) - t168;
t152 = t106 * t95;
t14 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t152 + t118;
t122 = qJDD(4) * qJ(5) + qJD(4) * t167 + t110 * t161 + t23;
t17 = -t114 * pkin(4) + t122;
t6 = t106 * t18 + t110 * t17;
t164 = -pkin(3) * t14 + pkin(8) * t6;
t79 = qJDD(4) + t87;
t156 = t106 * t79;
t85 = -t114 - t157;
t53 = t110 * t85 - t156;
t69 = -0.2e1 * t140 + t148;
t160 = pkin(3) * t69 + pkin(8) * t53;
t144 = t103 + t104;
t71 = t144 * t94;
t74 = t144 * t93;
t159 = pkin(3) * t74 + pkin(8) * t71;
t150 = t110 * t68;
t4 = t107 * t8 - t111 * t28;
t143 = pkin(2) * t4 + t165;
t142 = t106 * t28 - t170;
t141 = -t110 * t28 + t160;
t135 = t106 * (qJ(5) * t74 + t18) + t110 * ((-t114 + t74) * pkin(4) + t122) + t159;
t134 = t159 + t8;
t133 = t142 - t173;
t36 = t107 * t53 + t111 * t69;
t33 = pkin(2) * t36;
t132 = t33 + t141;
t126 = t107 * t93 - t111 * t94;
t131 = -pkin(2) * t126 - t136;
t42 = t107 * t71 + t111 * t74;
t41 = pkin(2) * t42;
t130 = t41 + t135;
t129 = t41 + t134;
t39 = t106 * t69 + t150;
t72 = -t107 * t94 - t111 * t93;
t116 = t152 * t167 - t118;
t123 = pkin(4) * t150 + t106 * (-pkin(4) * t140 + t116 + t168) + t170;
t121 = t69 * t145 + t110 * ((t69 - t140) * pkin(4) + t116) + t160;
t120 = t123 + t173;
t119 = t121 + t33;
t117 = t127 * t14 + t164;
t115 = pkin(2) * t72 - t32;
t75 = (t103 - t104) * t93;
t52 = t110 * (t114 - t158) + t156;
t46 = t69 * t110;
t45 = t68 * t106;
t30 = pkin(1) * (t108 * (-t107 * t74 + t111 * t71) + t112 * t42);
t19 = pkin(1) * (t108 * (-t107 * t69 + t111 * t53) + t112 * t36);
t12 = t107 * t32 - t111 * t136;
t11 = pkin(2) * t12;
t2 = t107 * t6 - t111 * t14;
t1 = pkin(2) * t2;
t3 = [0, 0, 0, 0, 0, qJDD(1), t138, t128, 0, 0, 0, 0, 0, 0, 0, t99, pkin(1) * (-t108 * t98 + t112 * t99) + t50, pkin(1) * (-t108 * t99 - t112 * t98) - t51, 0, pkin(1) * (t108 * t51 + t112 * t50), 0, 0, 0, 0, 0, t94, pkin(1) * (t108 * t72 - t112 * t126) + t131, pkin(1) * (t108 * t126 + t112 * t72) + t115, 0, pkin(1) * (t108 * (t107 * t136 + t111 * t32) + t112 * t12) + t11, t45, t39, t52, t46, t55, 0, t19 + t132, t133 - t174, t30 + t129, pkin(1) * (t108 * (t107 * t28 + t111 * t8) + t112 * t4) + t143, t45, t52, -t39, 0, -t55, t46, t119 + t19, t30 + t130, t120 + t174, pkin(1) * (t108 * (t107 * t14 + t111 * t6) + t112 * t2) + t1 - t14 * t163 - t14 * t145 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t50, -t51, 0, 0, 0, 0, 0, 0, 0, t94, t131, t115, 0, t11, t45, t39, t52, t46, t55, 0, t132, t133, t129, t143, t45, t52, -t39, 0, -t55, t46, t119, t130, t120, t1 + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t136, -t32, 0, 0, t45, t39, t52, t46, t55, 0, t141, t142, t134, t165, t45, t52, -t39, 0, -t55, t46, t121, t135, t123, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t75, t153, t87, t148, qJDD(4), -t22, -t23, 0, 0, -t87, t153, -t75, qJDD(4), -t148, t87, pkin(4) * t79 + qJ(5) * t85 - t18, (-pkin(4) * t106 + qJ(5) * t110) * t94, qJ(5) * t80 + (-t114 + t83) * pkin(4) + t122, -pkin(4) * t18 + qJ(5) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t153, -t83, t18;];
tauJ_reg = t3;
