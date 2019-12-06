% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:58
% EndTime: 2019-12-05 18:00:05
% DurationCPUTime: 1.39s
% Computational Cost: add. (3585->212), mult. (7380->256), div. (0->0), fcn. (4639->6), ass. (0->147)
t114 = qJDD(3) + qJDD(4);
t119 = sin(qJ(4));
t120 = sin(qJ(3));
t122 = cos(qJ(4));
t123 = cos(qJ(3));
t96 = (-t123 * t119 - t120 * t122) * qJD(1);
t157 = qJD(1) * t123;
t98 = -t119 * t120 * qJD(1) + t122 * t157;
t179 = t98 * t96;
t187 = t114 + t179;
t191 = t187 * pkin(4);
t190 = t119 * t187;
t189 = t122 * t187;
t188 = 2 * qJD(5);
t183 = pkin(6) + pkin(1);
t115 = qJD(3) + qJD(4);
t178 = t115 * t96;
t110 = t123 * qJDD(1);
t155 = qJD(1) * qJD(3);
t147 = t120 * t155;
t103 = t110 - t147;
t126 = qJD(1) ^ 2;
t159 = t123 * t126;
t121 = sin(qJ(1));
t124 = cos(qJ(1));
t145 = g(1) * t121 - t124 * g(2);
t137 = qJDD(2) - t145;
t171 = qJ(2) * t126;
t131 = t137 - t171;
t87 = -t183 * qJDD(1) + t131;
t172 = t123 * t87;
t57 = qJDD(3) * pkin(3) - pkin(7) * t103 + t172 + (-pkin(3) * t159 - pkin(7) * t155 + g(3)) * t120;
t146 = t123 * t155;
t153 = t120 * qJDD(1);
t102 = -t146 - t153;
t134 = qJD(3) * pkin(3) - pkin(7) * t157;
t117 = t120 ^ 2;
t163 = t117 * t126;
t75 = t123 * g(3) - t120 * t87;
t58 = -pkin(3) * t163 + t102 * pkin(7) - qJD(3) * t134 - t75;
t32 = t119 * t58 - t122 * t57;
t186 = -qJ(5) * t178 + t98 * t188 - t191 + t32;
t33 = t119 * t57 + t122 * t58;
t144 = -t122 * t102 + t103 * t119;
t62 = -qJD(4) * t98 - t144;
t82 = pkin(4) * t115 - qJ(5) * t98;
t94 = t96 ^ 2;
t16 = -t94 * pkin(4) + t62 * qJ(5) - t115 * t82 + t96 * t188 + t33;
t185 = -t98 * t82 - qJDD(5);
t116 = qJDD(1) * qJ(2);
t138 = g(1) * t124 + g(2) * t121;
t135 = -t116 + t138;
t184 = pkin(3) * t102 + (pkin(7) * t117 + t183) * t126 - t134 * t157 + t135;
t95 = t98 ^ 2;
t113 = t115 ^ 2;
t130 = (-qJD(4) + t115) * t98 - t144;
t136 = t102 * t119 + t103 * t122;
t63 = qJD(4) * t96 + t136;
t50 = t63 - t178;
t26 = t119 * t130 - t122 * t50;
t182 = pkin(7) * t26;
t67 = -t113 - t94;
t38 = t119 * t67 + t189;
t181 = pkin(7) * t38;
t69 = -t179 + t114;
t176 = t119 * t69;
t80 = -t95 - t113;
t51 = t122 * t80 - t176;
t180 = pkin(7) * t51;
t154 = qJD(2) * qJD(1);
t152 = -0.2e1 * t154;
t59 = t152 + t184;
t177 = t119 * t59;
t175 = t122 * t59;
t174 = t122 * t69;
t12 = t119 * t33 - t122 * t32;
t173 = t123 * t12;
t170 = qJ(5) * t119;
t169 = qJ(5) * t122;
t166 = qJDD(1) * pkin(1);
t165 = t115 * t119;
t164 = t115 * t122;
t118 = t123 ^ 2;
t162 = t118 * t126;
t151 = t120 * t159;
t161 = t120 * (qJDD(3) + t151);
t160 = t123 * (qJDD(3) - t151);
t158 = t117 + t118;
t156 = qJD(4) + t115;
t27 = t119 * t50 + t122 * t130;
t64 = -t94 - t95;
t150 = -pkin(3) * t64 + pkin(7) * t27;
t39 = t122 * t67 - t190;
t45 = t156 * t98 + t144;
t149 = -pkin(3) * t45 + pkin(7) * t39;
t48 = t156 * t96 + t136;
t52 = -t119 * t80 - t174;
t148 = -pkin(3) * t48 + pkin(7) * t52;
t13 = t119 * t32 + t122 * t33;
t9 = t120 * t27 + t123 * t26;
t143 = qJ(2) * t64 - t183 * t9;
t19 = t120 * t39 + t123 * t38;
t141 = qJ(2) * t45 - t183 * t19;
t28 = t120 * t52 + t123 * t51;
t140 = qJ(2) * t48 - t183 * t28;
t74 = g(3) * t120 + t172;
t53 = -t120 * t75 + t123 * t74;
t133 = pkin(4) * t80 - t16;
t15 = -qJ(5) * t63 - t186;
t128 = t15 + t191;
t111 = 0.2e1 * t154;
t127 = -pkin(4) * t62 + t111 - t184 - t185;
t125 = qJD(3) ^ 2;
t105 = t158 * qJDD(1);
t104 = t110 - 0.2e1 * t147;
t101 = 0.2e1 * t146 + t153;
t92 = -t131 + t166;
t85 = -t95 + t113;
t84 = t94 - t113;
t83 = t183 * t126 + t135 + t152;
t79 = -t161 + t123 * (-t125 - t162);
t78 = t120 * (-t125 - t163) + t160;
t71 = t95 - t94;
t49 = t63 + t178;
t44 = pkin(3) * t51;
t42 = pkin(4) * t50;
t37 = pkin(3) * t38;
t35 = (t123 * (t119 * t98 + t122 * t96) - t120 * (t119 * t96 - t122 * t98)) * t115;
t34 = -pkin(4) * t48 - qJ(5) * t69;
t30 = t123 * (t122 * t84 - t176) - t120 * (t119 * t84 + t174);
t29 = t123 * (-t119 * t85 + t189) - t120 * (t122 * t85 + t190);
t24 = pkin(3) * t26;
t22 = t123 * (t122 * t63 - t98 * t165) - t120 * (t119 * t63 + t98 * t164);
t21 = t123 * (-t119 * t62 - t96 * t164) - t120 * (t122 * t62 - t96 * t165);
t20 = -qJ(5) * t94 + t127;
t17 = (-t80 - t94) * qJ(5) + t127;
t14 = pkin(4) * t15;
t11 = (t67 + t94) * qJ(5) + (-t45 + t62) * pkin(4) + t59 + t185;
t10 = t123 * (-t119 * t49 - t122 * t45) - t120 * (-t119 * t45 + t122 * t49);
t7 = (t50 + t63) * qJ(5) + t186;
t6 = -pkin(4) * t64 + qJ(5) * t130 + t16;
t5 = -pkin(4) * t20 + qJ(5) * t16;
t4 = -t119 * t15 + t122 * t16;
t3 = t119 * t16 + t122 * t15;
t2 = t120 * t13 + t173;
t1 = t120 * t4 + t123 * t3;
t8 = [0, 0, 0, 0, 0, qJDD(1), t145, t138, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t137 - 0.2e1 * t166, t111 + 0.2e1 * t116 - t138, pkin(1) * t92 + qJ(2) * (-pkin(1) * t126 + t111 - t135), (t103 - t147) * t123, -t101 * t123 - t104 * t120, t160 - t120 * (t125 - t162), (-t102 + t146) * t120, t123 * (-t125 + t163) - t161, 0, qJ(2) * t101 - t120 * t83 - t183 * t78, qJ(2) * t104 - t123 * t83 - t183 * t79, t183 * t105 - t158 * t171 - t53, -qJ(2) * t83 - t183 * t53, t22, t10, t29, t21, t30, t35, t123 * (-t177 - t181) - t120 * (t149 + t175) + t141, t123 * (-t175 - t180) - t120 * (t148 - t177) + t140, t123 * (-t12 - t182) - t120 * (t13 + t150) + t143, -pkin(7) * t173 - t120 * (pkin(3) * t59 + pkin(7) * t13) - qJ(2) * t59 - t183 * t2, t22, t10, t29, t21, t30, t35, t123 * (-t11 * t119 - t169 * t187 - t181) - t120 * (t11 * t122 - t170 * t187 + t149) + t141, t123 * (-t119 * t34 + t122 * t17 - t180) - t120 * (t119 * t17 + t122 * t34 + t148) + t140, t123 * (-t119 * t6 + t122 * t7 - t182) - t120 * (t119 * t7 + t122 * t6 + t150) + t143, t123 * (-pkin(7) * t3 - t119 * t5 - t15 * t169) - t120 * (-pkin(3) * t20 + pkin(7) * t4 + t122 * t5 - t15 * t170) + qJ(2) * t20 - t183 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t126, -t92, 0, 0, 0, 0, 0, 0, t78, t79, -t105, t53, 0, 0, 0, 0, 0, 0, t19, t28, t9, t2, 0, 0, 0, 0, 0, 0, t19, t28, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, (-t117 + t118) * t126, t110, -t151, -t153, qJDD(3), t74, t75, 0, 0, -t179, t71, t50, t179, t130, t114, -t32 + t37, t44 - t33, t24, pkin(3) * t12, -t179, t71, t50, t179, t130, t114, t128 + t37, t44 + t133, -t42 + t24, pkin(3) * t3 + t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, t71, t50, t179, t130, t114, -t32, -t33, 0, 0, -t179, t71, t50, t179, t130, t114, t128, t133, -t42, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t48, t64, t20;];
tauJ_reg = t8;
