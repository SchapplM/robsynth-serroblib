% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:41
% EndTime: 2019-12-31 20:53:45
% DurationCPUTime: 0.90s
% Computational Cost: add. (2482->177), mult. (3433->179), div. (0->0), fcn. (1597->6), ass. (0->111)
t117 = sin(qJ(2));
t120 = cos(qJ(2));
t116 = sin(qJ(3));
t119 = cos(qJ(3));
t109 = qJD(1) + qJD(2);
t107 = t109 ^ 2;
t149 = t119 * t107 * t116;
t85 = qJDD(3) - t149;
t164 = t119 * t85;
t122 = qJD(3) ^ 2;
t113 = t116 ^ 2;
t159 = t113 * t107;
t89 = -t122 - t159;
t56 = t116 * t89 + t164;
t153 = qJD(3) * t109;
t146 = t119 * t153;
t108 = qJDD(1) + qJDD(2);
t156 = t116 * t108;
t69 = 0.2e1 * t146 + t156;
t28 = pkin(1) * (t117 * t56 + t120 * t69);
t84 = qJDD(3) + t149;
t168 = t116 * t84;
t114 = t119 ^ 2;
t158 = t114 * t107;
t91 = -t122 - t158;
t55 = -t119 * t91 + t168;
t155 = t119 * t108;
t99 = t116 * t153;
t72 = -0.2e1 * t99 + t155;
t27 = pkin(1) * (t117 * t55 - t120 * t72);
t173 = pkin(2) * t69 + pkin(7) * t56;
t70 = t146 + t156;
t186 = t70 + t146;
t177 = t119 * g(3);
t132 = -qJDD(3) * pkin(3) - t122 * qJ(4) + qJDD(4) + t177;
t124 = -(2 * qJD(5) * qJD(3)) + t132 - qJ(5) * t84 + (-t146 + t70) * pkin(4);
t157 = t116 * qJ(4);
t134 = -t119 * pkin(3) - t157;
t169 = t107 * t134;
t118 = sin(qJ(1));
t121 = cos(qJ(1));
t142 = t118 * g(1) - t121 * g(2);
t130 = qJDD(1) * pkin(1) + t142;
t138 = t121 * g(1) + t118 * g(2);
t82 = -qJD(1) ^ 2 * pkin(1) - t138;
t43 = t117 * t130 + t120 * t82;
t37 = -t107 * pkin(2) + t108 * pkin(7) + t43;
t141 = t37 + t169;
t135 = t141 * t116;
t13 = t135 + t124;
t174 = pkin(2) * t72 - pkin(7) * t55;
t184 = (-t72 + t99) * pkin(3);
t30 = t116 * t37 + t177;
t31 = -t116 * g(3) + t119 * t37;
t18 = t116 * t30 + t119 * t31;
t183 = -pkin(3) * t89 + qJ(4) * t85;
t160 = t109 * t116;
t83 = pkin(4) * t160 - qJD(3) * qJ(5);
t182 = -pkin(4) * t158 - t83 * t160;
t50 = t164 + t116 * (-t122 + t158);
t181 = qJD(3) * t83 + qJDD(5);
t42 = -t117 * t82 + t120 * t130;
t36 = -t108 * pkin(2) - t107 * pkin(7) - t42;
t71 = -t99 + t155;
t128 = -t71 * pkin(3) - t186 * qJ(4) + t36;
t133 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t116;
t20 = t109 * t133 + t128;
t136 = qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t119 * t169 + t31;
t23 = -t122 * pkin(3) + t136;
t24 = t135 + t132;
t8 = t116 * t24 + t119 * t23;
t178 = -pkin(2) * t20 + pkin(7) * t8;
t176 = pkin(3) + qJ(5);
t175 = -pkin(2) * t36 + pkin(7) * t18;
t154 = t113 + t114;
t75 = t154 * t108;
t78 = t154 * t107;
t172 = pkin(2) * t78 + pkin(7) * t75;
t171 = qJ(4) * t78;
t170 = qJ(4) * t91;
t165 = t119 * t69;
t161 = t71 * qJ(5);
t152 = qJD(5) * t119;
t151 = t116 * t36 - t173;
t150 = -t119 * t36 + t174;
t126 = t128 + t182;
t11 = -t161 + (t133 - 0.2e1 * t152) * t109 + t126;
t125 = t71 * pkin(4) + t181 + t23;
t16 = -qJ(5) * t158 + t125;
t4 = t116 * t13 + t119 * t16;
t148 = t116 * (pkin(4) * t13 - qJ(4) * t11) + t119 * (pkin(4) * t16 - t176 * t11) - pkin(2) * t11 + pkin(7) * t4;
t147 = qJD(4) * t160;
t127 = (-t122 + t78) * pkin(3) + t136;
t145 = t116 * (t171 + (pkin(4) * t108 + t141) * t116 + t124) + t119 * ((t78 - t158) * qJ(5) + (t71 + t155) * pkin(4) + t127 + t181) + t172;
t96 = 0.2e1 * t109 * t152;
t97 = 0.2e1 * t147;
t144 = t116 * (-pkin(4) * t84 + qJ(4) * t72) + t119 * (pkin(4) * t91 + t96 + t97 + (t71 + t72) * qJ(5) - t184 - t126) + t174;
t123 = -pkin(3) * t99 + qJ(4) * t69 - t128 + t97;
t143 = t116 * (pkin(4) * t89 + t123 + t161 - t182 + t96) + t119 * (pkin(4) * t85 + t176 * t69) + t173;
t140 = t119 * t127 + t172 + t116 * (t24 + t171);
t139 = t172 + t18;
t38 = t116 * t72 + t165;
t52 = t119 * (t122 - t159) + t168;
t131 = pkin(3) * t165 + t116 * t123 + t173;
t129 = -t72 * t157 + t119 * (t128 - 0.2e1 * t147 + t184) - t174;
t98 = qJ(4) * t155;
t79 = (t113 - t114) * t107;
t45 = t186 * t116;
t44 = (t71 - t99) * t119;
t40 = pkin(1) * (t117 * t75 + t120 * t78);
t1 = [0, 0, 0, 0, 0, qJDD(1), t142, t138, 0, 0, 0, 0, 0, 0, 0, t108, pkin(1) * (-t117 * t107 + t120 * t108) + t42, pkin(1) * (-t120 * t107 - t117 * t108) - t43, 0, pkin(1) * (t117 * t43 + t120 * t42), t45, t38, t52, t44, t50, 0, -t27 + t150, t151 - t28, t40 + t139, pkin(1) * (t117 * t18 - t120 * t36) + t175, 0, -t52, -t50, t45, t38, t44, t40 + t140, t129 + t27, t131 + t28, pkin(1) * t117 * t8 + (-pkin(1) * t120 + t134) * t20 + t178, 0, -t50, t52, t44, -t38, t45, t40 + t145, t28 + t143, -t27 + t144, pkin(1) * (-t120 * t11 + t117 * t4) + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t42, -t43, 0, 0, t45, t38, t52, t44, t50, 0, t150, t151, t139, t175, 0, -t52, -t50, t45, t38, t44, t140, t129, t131, t134 * t20 + t178, 0, -t50, t52, t44, -t38, t45, t145, t143, t144, t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t79, t156, t149, t155, qJDD(3), -t30, -t31, 0, 0, qJDD(3), -t156, -t155, -t149, t79, t149, -pkin(3) * t156 + t98, -pkin(3) * t84 - t170 + t24, t23 + t183, -pkin(3) * t24 + qJ(4) * t23, qJDD(3), -t155, t156, t149, -t79, -t149, -t176 * t156 + t98, (-t89 - t158) * qJ(5) + t125 + t183, t176 * t84 - t13 + t170, qJ(4) * t16 - t176 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t84, t89, t24, 0, 0, 0, 0, 0, 0, t156, t89, -t84, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t85, t91, t16;];
tauJ_reg = t1;
