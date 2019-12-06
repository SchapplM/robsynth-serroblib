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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:49:39
% EndTime: 2019-12-05 17:49:44
% DurationCPUTime: 1.19s
% Computational Cost: add. (4646->177), mult. (6931->273), div. (0->0), fcn. (4288->10), ass. (0->121)
t116 = qJDD(1) + qJDD(3);
t155 = (qJD(1) + qJD(3));
t153 = t155 ^ 2;
t123 = sin(qJ(3));
t126 = cos(qJ(3));
t119 = sin(pkin(8));
t121 = cos(pkin(8));
t129 = qJD(1) ^ 2;
t124 = sin(qJ(1));
t127 = cos(qJ(1));
t142 = g(2) * t124 - g(3) * t127;
t138 = -pkin(1) * t129 + t142;
t158 = t127 * g(2) + t124 * g(3);
t166 = qJDD(1) * pkin(1);
t141 = t158 + t166;
t135 = -t119 * t138 + t121 * t141;
t134 = qJDD(1) * pkin(2) + t135;
t171 = t119 * t141 + t121 * t138;
t73 = -pkin(2) * t129 + t171;
t54 = t123 * t134 + t126 * t73;
t178 = -t153 * pkin(3) + t116 * qJ(4) + (2 * qJD(4) * t155) + t54;
t122 = sin(qJ(5));
t118 = sin(pkin(9));
t120 = cos(pkin(9));
t125 = cos(qJ(5));
t147 = t125 * t155;
t148 = t122 * t155;
t84 = t118 * t148 - t120 * t147;
t86 = t118 * t147 + t120 * t148;
t70 = t86 * t84;
t174 = qJDD(5) - t70;
t177 = t122 * t174;
t176 = t125 * t174;
t117 = -g(1) + qJDD(2);
t111 = t120 * t117;
t41 = t178 * t118 - t111;
t42 = t118 * t117 + t178 * t120;
t18 = t118 * t41 + t120 * t42;
t175 = t111 + (t153 * t120 * pkin(4) - t116 * pkin(7) - t178) * t118;
t130 = t118 ^ 2;
t132 = t120 ^ 2;
t145 = t132 * t153;
t96 = t130 * t153 + t145;
t82 = t84 ^ 2;
t83 = t86 ^ 2;
t164 = t120 * t116;
t35 = -pkin(4) * t145 + pkin(7) * t164 + t42;
t13 = t122 * t35 - t125 * t175;
t14 = t122 * t175 + t125 * t35;
t7 = t122 * t14 - t125 * t13;
t173 = t118 * t7;
t53 = -t123 * t73 + t126 * t134;
t51 = -t116 * pkin(3) - t153 * qJ(4) + qJDD(4) - t53;
t172 = -pkin(3) * t51 + qJ(4) * t18;
t37 = -pkin(4) * t164 - pkin(7) * t96 + t51;
t170 = t122 * t37;
t64 = qJDD(5) + t70;
t169 = t122 * t64;
t168 = t125 * t37;
t167 = t125 * t64;
t165 = t118 * t116;
t163 = t123 * t116;
t162 = t126 * t116;
t161 = t84 * qJD(5);
t160 = t86 * qJD(5);
t157 = qJD(5) * t122;
t156 = qJD(5) * t125;
t90 = t96 * t120;
t154 = pkin(3) * t164 - qJ(4) * t90 - t120 * t51;
t57 = -t122 * t165 + t125 * t164;
t81 = (t118 * t125 + t120 * t122) * t116;
t47 = t122 * t57 - t125 * t81;
t48 = t122 * t81 + t125 * t57;
t26 = -t118 * t47 + t120 * t48;
t59 = -t82 - t83;
t8 = t122 * t13 + t125 * t14;
t152 = t118 * (-pkin(7) * t47 - t7) + t120 * (-pkin(4) * t59 + pkin(7) * t48 + t8) - pkin(3) * t59 + qJ(4) * t26;
t128 = qJD(5) ^ 2;
t62 = -t128 - t82;
t43 = t122 * t62 + t176;
t44 = t125 * t62 - t177;
t24 = -t118 * t43 + t120 * t44;
t66 = -t57 + 0.2e1 * t160;
t151 = t118 * (-pkin(7) * t43 + t170) + t120 * (-pkin(4) * t66 + pkin(7) * t44 - t168) - pkin(3) * t66 + qJ(4) * t24;
t78 = -t83 - t128;
t55 = t125 * t78 - t169;
t56 = -t122 * t78 - t167;
t31 = -t118 * t55 + t120 * t56;
t68 = t81 - 0.2e1 * t161;
t150 = t118 * (-pkin(7) * t55 + t168) + t120 * (-pkin(4) * t68 + pkin(7) * t56 + t170) - pkin(3) * t68 + qJ(4) * t31;
t108 = t130 * t116;
t109 = t132 * t116;
t94 = t109 + t108;
t149 = pkin(3) * t96 + qJ(4) * t94 + t18;
t89 = t96 * t118;
t140 = -pkin(3) * t165 + qJ(4) * t89 + t118 * t51;
t3 = t120 * t8 - t173;
t139 = -pkin(7) * t173 + qJ(4) * t3 - pkin(3) * t37 + t120 * (-pkin(4) * t37 + pkin(7) * t8);
t137 = t123 * t153 - t162;
t97 = -t126 * t153 - t163;
t99 = 0.2e1 * t118 * t164;
t77 = -t83 + t128;
t76 = t82 - t128;
t75 = t120 * t162 - t123 * t90;
t74 = -t118 * t162 + t123 * t89;
t72 = t123 * t94 + t126 * t96;
t69 = t81 - t161;
t67 = t57 - t160;
t40 = (t118 * (t122 * t86 - t125 * t84) + t120 * (-t122 * t84 - t125 * t86)) * qJD(5);
t33 = t118 * (t125 * t69 - t86 * t157) + t120 * (t122 * t69 + t86 * t156);
t32 = t118 * (-t122 * t67 + t84 * t156) + t120 * (t125 * t67 + t84 * t157);
t30 = t118 * (-t122 * t77 + t176) + t120 * (t125 * t77 + t177);
t29 = t118 * (t125 * t76 - t169) + t120 * (t122 * t76 + t167);
t27 = t123 * t54 + t126 * t53;
t25 = t118 * (-t122 * t68 - t125 * t66) + t120 * (-t122 * t66 + t125 * t68);
t20 = t123 * t31 - t126 * t68;
t15 = t123 * t24 - t126 * t66;
t10 = t123 * t26 - t126 * t59;
t9 = t123 * t18 - t126 * t51;
t1 = t123 * t3 - t126 * t37;
t2 = [0, 0, 0, 0, 0, qJDD(1), t158, -t142, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t119 * t142 + (t158 + 0.2e1 * t166) * t121, pkin(1) * (-qJDD(1) * t119 - t121 * t129) - t171, 0, pkin(1) * (t119 * t171 + t121 * t135), 0, 0, 0, 0, 0, t116, pkin(1) * (t119 * t97 - t121 * t137) - pkin(2) * t137 + t53, pkin(1) * (t119 * t137 + t121 * t97) + pkin(2) * t97 - t54, 0, pkin(1) * (t119 * (-t123 * t53 + t126 * t54) + t121 * t27) + pkin(2) * t27, t108, t99, 0, t109, 0, 0, pkin(1) * (t119 * (-t120 * t163 - t126 * t90) + t121 * t75) + pkin(2) * t75 + t154, pkin(1) * (t119 * (t118 * t163 + t126 * t89) + t121 * t74) + pkin(2) * t74 + t140, pkin(1) * (t119 * (-t123 * t96 + t126 * t94) + t121 * t72) + pkin(2) * t72 + t149, pkin(1) * (t119 * (t123 * t51 + t126 * t18) + t121 * t9) + pkin(2) * t9 + t172, t33, t25, t30, t32, t29, t40, pkin(1) * (t119 * (t123 * t66 + t126 * t24) + t121 * t15) + pkin(2) * t15 + t151, pkin(1) * (t119 * (t123 * t68 + t126 * t31) + t121 * t20) + pkin(2) * t20 + t150, pkin(1) * (t119 * (t123 * t59 + t126 * t26) + t121 * t10) + pkin(2) * t10 + t152, pkin(1) * (t119 * (t123 * t37 + t126 * t3) + t121 * t1) + pkin(2) * t1 + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118 * t42 - t120 * t41, 0, 0, 0, 0, 0, 0, t118 * t44 + t120 * t43, t118 * t56 + t120 * t55, t118 * t48 + t120 * t47, t118 * t8 + t120 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t53, -t54, 0, 0, t108, t99, 0, t109, 0, 0, t154, t140, t149, t172, t33, t25, t30, t32, t29, t40, t151, t150, t152, t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t165, -t96, t51, 0, 0, 0, 0, 0, 0, t66, t68, t59, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t83 - t82, t81, -t70, t57, qJDD(5), -t13, -t14, 0, 0;];
tauJ_reg = t2;
