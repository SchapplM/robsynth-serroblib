% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:07
% EndTime: 2022-01-23 09:19:10
% DurationCPUTime: 1.17s
% Computational Cost: add. (4646->176), mult. (6931->273), div. (0->0), fcn. (4288->10), ass. (0->121)
t115 = qJDD(1) + qJDD(3);
t155 = (qJD(1) + qJD(3));
t153 = t155 ^ 2;
t122 = sin(qJ(3));
t125 = cos(qJ(3));
t118 = sin(pkin(8));
t120 = cos(pkin(8));
t128 = qJD(1) ^ 2;
t123 = sin(qJ(1));
t126 = cos(qJ(1));
t141 = g(1) * t126 + g(2) * t123;
t137 = -pkin(1) * t128 - t141;
t151 = t123 * g(1) - g(2) * t126;
t165 = qJDD(1) * pkin(1);
t138 = t151 + t165;
t134 = -t118 * t137 + t120 * t138;
t133 = qJDD(1) * pkin(2) + t134;
t170 = t118 * t138 + t120 * t137;
t73 = -pkin(2) * t128 + t170;
t54 = t122 * t133 + t125 * t73;
t177 = -t153 * pkin(3) + t115 * qJ(4) + (2 * qJD(4) * t155) + t54;
t121 = sin(qJ(5));
t117 = sin(pkin(9));
t119 = cos(pkin(9));
t124 = cos(qJ(5));
t146 = t124 * t155;
t147 = t121 * t155;
t84 = t117 * t147 - t119 * t146;
t86 = t117 * t146 + t119 * t147;
t70 = t86 * t84;
t173 = qJDD(5) - t70;
t176 = t121 * t173;
t175 = t124 * t173;
t116 = -g(3) + qJDD(2);
t111 = t119 * t116;
t41 = t177 * t117 - t111;
t42 = t117 * t116 + t177 * t119;
t18 = t117 * t41 + t119 * t42;
t174 = t111 + (t153 * t119 * pkin(4) - t115 * pkin(7) - t177) * t117;
t129 = t117 ^ 2;
t131 = t119 ^ 2;
t144 = t131 * t153;
t96 = t129 * t153 + t144;
t82 = t84 ^ 2;
t83 = t86 ^ 2;
t163 = t119 * t115;
t35 = -pkin(4) * t144 + pkin(7) * t163 + t42;
t13 = t121 * t35 - t124 * t174;
t14 = t174 * t121 + t124 * t35;
t7 = t121 * t14 - t124 * t13;
t172 = t117 * t7;
t53 = -t122 * t73 + t125 * t133;
t51 = -t115 * pkin(3) - t153 * qJ(4) + qJDD(4) - t53;
t171 = -pkin(3) * t51 + qJ(4) * t18;
t37 = -pkin(4) * t163 - t96 * pkin(7) + t51;
t169 = t121 * t37;
t64 = qJDD(5) + t70;
t168 = t121 * t64;
t167 = t124 * t37;
t166 = t124 * t64;
t164 = t117 * t115;
t162 = t122 * t115;
t161 = t125 * t115;
t160 = t84 * qJD(5);
t159 = t86 * qJD(5);
t157 = qJD(5) * t121;
t156 = qJD(5) * t124;
t90 = t96 * t119;
t154 = pkin(3) * t163 - qJ(4) * t90 - t119 * t51;
t57 = -t121 * t164 + t124 * t163;
t81 = (t124 * t117 + t121 * t119) * t115;
t47 = t121 * t57 - t124 * t81;
t48 = t121 * t81 + t124 * t57;
t26 = -t117 * t47 + t119 * t48;
t59 = -t82 - t83;
t8 = t121 * t13 + t124 * t14;
t152 = t117 * (-pkin(7) * t47 - t7) + t119 * (-pkin(4) * t59 + pkin(7) * t48 + t8) - pkin(3) * t59 + qJ(4) * t26;
t127 = qJD(5) ^ 2;
t62 = -t127 - t82;
t43 = t121 * t62 + t175;
t44 = t124 * t62 - t176;
t24 = -t117 * t43 + t119 * t44;
t66 = -t57 + 0.2e1 * t159;
t150 = t117 * (-pkin(7) * t43 + t169) + t119 * (-pkin(4) * t66 + pkin(7) * t44 - t167) - pkin(3) * t66 + qJ(4) * t24;
t78 = -t83 - t127;
t55 = t124 * t78 - t168;
t56 = -t121 * t78 - t166;
t31 = -t117 * t55 + t119 * t56;
t68 = t81 - 0.2e1 * t160;
t149 = t117 * (-pkin(7) * t55 + t167) + t119 * (-pkin(4) * t68 + pkin(7) * t56 + t169) - pkin(3) * t68 + qJ(4) * t31;
t108 = t129 * t115;
t109 = t131 * t115;
t94 = t109 + t108;
t148 = pkin(3) * t96 + qJ(4) * t94 + t18;
t89 = t96 * t117;
t140 = -pkin(3) * t164 + qJ(4) * t89 + t117 * t51;
t3 = t119 * t8 - t172;
t139 = -pkin(7) * t172 + qJ(4) * t3 - pkin(3) * t37 + t119 * (-pkin(4) * t37 + pkin(7) * t8);
t136 = t122 * t153 - t161;
t97 = -t125 * t153 - t162;
t99 = 0.2e1 * t117 * t163;
t77 = -t83 + t127;
t76 = t82 - t127;
t75 = t119 * t161 - t122 * t90;
t74 = -t117 * t161 + t122 * t89;
t72 = t122 * t94 + t125 * t96;
t69 = t81 - t160;
t67 = t57 - t159;
t40 = (t117 * (t121 * t86 - t124 * t84) + t119 * (-t121 * t84 - t124 * t86)) * qJD(5);
t33 = t117 * (t124 * t69 - t86 * t157) + t119 * (t121 * t69 + t86 * t156);
t32 = t117 * (-t121 * t67 + t84 * t156) + t119 * (t124 * t67 + t84 * t157);
t30 = t117 * (-t121 * t77 + t175) + t119 * (t124 * t77 + t176);
t29 = t117 * (t124 * t76 - t168) + t119 * (t121 * t76 + t166);
t27 = t122 * t54 + t125 * t53;
t25 = t117 * (-t121 * t68 - t124 * t66) + t119 * (-t121 * t66 + t124 * t68);
t20 = t122 * t31 - t125 * t68;
t15 = t122 * t24 - t125 * t66;
t10 = t122 * t26 - t125 * t59;
t9 = t122 * t18 - t125 * t51;
t1 = t122 * t3 - t125 * t37;
t2 = [0, 0, 0, 0, 0, qJDD(1), t151, t141, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t118 * t141 + (t151 + 0.2e1 * t165) * t120, pkin(1) * (-qJDD(1) * t118 - t120 * t128) - t170, 0, pkin(1) * (t118 * t170 + t120 * t134), 0, 0, 0, 0, 0, t115, pkin(1) * (t118 * t97 - t120 * t136) - pkin(2) * t136 + t53, pkin(1) * (t118 * t136 + t120 * t97) + pkin(2) * t97 - t54, 0, pkin(1) * (t118 * (-t122 * t53 + t125 * t54) + t120 * t27) + pkin(2) * t27, t108, t99, 0, t109, 0, 0, pkin(1) * (t118 * (-t119 * t162 - t125 * t90) + t120 * t75) + pkin(2) * t75 + t154, pkin(1) * (t118 * (t117 * t162 + t125 * t89) + t120 * t74) + pkin(2) * t74 + t140, pkin(1) * (t118 * (-t122 * t96 + t125 * t94) + t120 * t72) + pkin(2) * t72 + t148, pkin(1) * (t118 * (t122 * t51 + t125 * t18) + t120 * t9) + pkin(2) * t9 + t171, t33, t25, t30, t32, t29, t40, pkin(1) * (t118 * (t122 * t66 + t125 * t24) + t120 * t15) + pkin(2) * t15 + t150, pkin(1) * (t118 * (t122 * t68 + t125 * t31) + t120 * t20) + pkin(2) * t20 + t149, pkin(1) * (t118 * (t122 * t59 + t125 * t26) + t120 * t10) + pkin(2) * t10 + t152, pkin(1) * (t118 * (t122 * t37 + t125 * t3) + t120 * t1) + pkin(2) * t1 + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117 * t42 - t119 * t41, 0, 0, 0, 0, 0, 0, t117 * t44 + t119 * t43, t117 * t56 + t119 * t55, t117 * t48 + t119 * t47, t117 * t8 + t119 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t53, -t54, 0, 0, t108, t99, 0, t109, 0, 0, t154, t140, t148, t171, t33, t25, t30, t32, t29, t40, t150, t149, t152, t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, t164, -t96, t51, 0, 0, 0, 0, 0, 0, t66, t68, t59, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t83 - t82, t81, -t70, t57, qJDD(5), -t13, -t14, 0, 0;];
tauJ_reg = t2;
