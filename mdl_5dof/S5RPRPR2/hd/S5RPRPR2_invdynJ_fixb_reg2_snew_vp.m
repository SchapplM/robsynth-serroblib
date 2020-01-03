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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:34:03
% EndTime: 2020-01-03 11:34:10
% DurationCPUTime: 1.25s
% Computational Cost: add. (4646->176), mult. (6931->273), div. (0->0), fcn. (4288->10), ass. (0->121)
t114 = qJDD(1) + qJDD(3);
t154 = (qJD(1) + qJD(3));
t152 = t154 ^ 2;
t121 = sin(qJ(3));
t124 = cos(qJ(3));
t117 = sin(pkin(8));
t119 = cos(pkin(8));
t122 = sin(qJ(1));
t125 = cos(qJ(1));
t141 = -g(2) * t125 - g(3) * t122;
t164 = qJDD(1) * pkin(1);
t135 = t141 + t164;
t127 = qJD(1) ^ 2;
t140 = g(2) * t122 - g(3) * t125;
t137 = -pkin(1) * t127 - t140;
t133 = -t117 * t137 + t119 * t135;
t132 = qJDD(1) * pkin(2) + t133;
t169 = t117 * t135 + t119 * t137;
t73 = -pkin(2) * t127 + t169;
t54 = t121 * t132 + t124 * t73;
t176 = -t152 * pkin(3) + t114 * qJ(4) + (2 * qJD(4) * t154) + t54;
t120 = sin(qJ(5));
t116 = sin(pkin(9));
t118 = cos(pkin(9));
t123 = cos(qJ(5));
t146 = t123 * t154;
t147 = t120 * t154;
t84 = t116 * t147 - t118 * t146;
t86 = t116 * t146 + t118 * t147;
t70 = t86 * t84;
t172 = qJDD(5) - t70;
t175 = t120 * t172;
t174 = t123 * t172;
t115 = -g(1) + qJDD(2);
t111 = t118 * t115;
t41 = t176 * t116 - t111;
t42 = t116 * t115 + t176 * t118;
t18 = t116 * t41 + t118 * t42;
t173 = t111 + (t152 * t118 * pkin(4) - t114 * pkin(7) - t176) * t116;
t128 = t116 ^ 2;
t130 = t118 ^ 2;
t144 = t130 * t152;
t96 = t128 * t152 + t144;
t82 = t84 ^ 2;
t83 = t86 ^ 2;
t162 = t118 * t114;
t35 = -pkin(4) * t144 + pkin(7) * t162 + t42;
t13 = t120 * t35 - t123 * t173;
t14 = t120 * t173 + t123 * t35;
t7 = t120 * t14 - t123 * t13;
t171 = t116 * t7;
t53 = -t121 * t73 + t124 * t132;
t51 = -t114 * pkin(3) - t152 * qJ(4) + qJDD(4) - t53;
t170 = -pkin(3) * t51 + qJ(4) * t18;
t37 = -pkin(4) * t162 - pkin(7) * t96 + t51;
t168 = t120 * t37;
t64 = qJDD(5) + t70;
t167 = t120 * t64;
t166 = t123 * t37;
t165 = t123 * t64;
t163 = t116 * t114;
t161 = t121 * t114;
t160 = t124 * t114;
t159 = t84 * qJD(5);
t158 = t86 * qJD(5);
t156 = qJD(5) * t120;
t155 = qJD(5) * t123;
t90 = t96 * t118;
t153 = pkin(3) * t162 - qJ(4) * t90 - t118 * t51;
t57 = -t120 * t163 + t123 * t162;
t81 = (t123 * t116 + t120 * t118) * t114;
t47 = t120 * t57 - t123 * t81;
t48 = t120 * t81 + t123 * t57;
t26 = -t116 * t47 + t118 * t48;
t59 = -t82 - t83;
t8 = t120 * t13 + t123 * t14;
t151 = t116 * (-pkin(7) * t47 - t7) + t118 * (-pkin(4) * t59 + pkin(7) * t48 + t8) - pkin(3) * t59 + qJ(4) * t26;
t126 = qJD(5) ^ 2;
t62 = -t126 - t82;
t43 = t120 * t62 + t174;
t44 = t123 * t62 - t175;
t24 = -t116 * t43 + t118 * t44;
t66 = -t57 + 0.2e1 * t158;
t150 = t116 * (-pkin(7) * t43 + t168) + t118 * (-pkin(4) * t66 + pkin(7) * t44 - t166) - pkin(3) * t66 + qJ(4) * t24;
t78 = -t83 - t126;
t55 = t123 * t78 - t167;
t56 = -t120 * t78 - t165;
t31 = -t116 * t55 + t118 * t56;
t68 = t81 - 0.2e1 * t159;
t149 = t116 * (-pkin(7) * t55 + t166) + t118 * (-pkin(4) * t68 + pkin(7) * t56 + t168) - pkin(3) * t68 + qJ(4) * t31;
t108 = t128 * t114;
t109 = t130 * t114;
t94 = t109 + t108;
t148 = pkin(3) * t96 + qJ(4) * t94 + t18;
t89 = t96 * t116;
t139 = -pkin(3) * t163 + qJ(4) * t89 + t116 * t51;
t3 = t118 * t8 - t171;
t138 = -pkin(7) * t171 + qJ(4) * t3 - pkin(3) * t37 + t118 * (-pkin(4) * t37 + pkin(7) * t8);
t136 = t121 * t152 - t160;
t97 = -t124 * t152 - t161;
t99 = 0.2e1 * t116 * t162;
t77 = -t83 + t126;
t76 = t82 - t126;
t75 = t118 * t160 - t121 * t90;
t74 = -t116 * t160 + t121 * t89;
t72 = t121 * t94 + t124 * t96;
t69 = t81 - t159;
t67 = t57 - t158;
t40 = (t116 * (t120 * t86 - t123 * t84) + t118 * (-t120 * t84 - t123 * t86)) * qJD(5);
t33 = t116 * (t123 * t69 - t86 * t156) + t118 * (t120 * t69 + t86 * t155);
t32 = t116 * (-t120 * t67 + t84 * t155) + t118 * (t123 * t67 + t84 * t156);
t30 = t116 * (-t120 * t77 + t174) + t118 * (t123 * t77 + t175);
t29 = t116 * (t123 * t76 - t167) + t118 * (t120 * t76 + t165);
t27 = t121 * t54 + t124 * t53;
t25 = t116 * (-t120 * t68 - t123 * t66) + t118 * (-t120 * t66 + t123 * t68);
t20 = t121 * t31 - t124 * t68;
t15 = t121 * t24 - t124 * t66;
t10 = t121 * t26 - t124 * t59;
t9 = t121 * t18 - t124 * t51;
t1 = t121 * t3 - t124 * t37;
t2 = [0, 0, 0, 0, 0, qJDD(1), t141, t140, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t117 * t140 + (t141 + 0.2e1 * t164) * t119, pkin(1) * (-qJDD(1) * t117 - t119 * t127) - t169, 0, pkin(1) * (t117 * t169 + t119 * t133), 0, 0, 0, 0, 0, t114, pkin(1) * (t117 * t97 - t119 * t136) - pkin(2) * t136 + t53, pkin(1) * (t117 * t136 + t119 * t97) + pkin(2) * t97 - t54, 0, pkin(1) * (t117 * (-t121 * t53 + t124 * t54) + t119 * t27) + pkin(2) * t27, t108, t99, 0, t109, 0, 0, pkin(1) * (t117 * (-t118 * t161 - t124 * t90) + t119 * t75) + pkin(2) * t75 + t153, pkin(1) * (t117 * (t116 * t161 + t124 * t89) + t119 * t74) + pkin(2) * t74 + t139, pkin(1) * (t117 * (-t121 * t96 + t124 * t94) + t119 * t72) + pkin(2) * t72 + t148, pkin(1) * (t117 * (t121 * t51 + t124 * t18) + t119 * t9) + pkin(2) * t9 + t170, t33, t25, t30, t32, t29, t40, pkin(1) * (t117 * (t121 * t66 + t124 * t24) + t119 * t15) + pkin(2) * t15 + t150, pkin(1) * (t117 * (t121 * t68 + t124 * t31) + t119 * t20) + pkin(2) * t20 + t149, pkin(1) * (t117 * (t121 * t59 + t124 * t26) + t119 * t10) + pkin(2) * t10 + t151, pkin(1) * (t117 * (t121 * t37 + t124 * t3) + t119 * t1) + pkin(2) * t1 + t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116 * t42 - t118 * t41, 0, 0, 0, 0, 0, 0, t116 * t44 + t118 * t43, t116 * t56 + t118 * t55, t116 * t48 + t118 * t47, t116 * t8 + t118 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t53, -t54, 0, 0, t108, t99, 0, t109, 0, 0, t153, t139, t148, t170, t33, t25, t30, t32, t29, t40, t150, t149, t151, t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, t163, -t96, t51, 0, 0, 0, 0, 0, 0, t66, t68, t59, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t83 - t82, t81, -t70, t57, qJDD(5), -t13, -t14, 0, 0;];
tauJ_reg = t2;
