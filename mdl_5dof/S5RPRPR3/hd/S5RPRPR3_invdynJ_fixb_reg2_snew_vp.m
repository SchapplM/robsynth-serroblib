% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR3
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:33
% EndTime: 2019-12-05 17:51:38
% DurationCPUTime: 1.03s
% Computational Cost: add. (4082->165), mult. (6022->245), div. (0->0), fcn. (3535->10), ass. (0->123)
t117 = sin(pkin(9));
t114 = qJD(1) + qJD(3);
t121 = sin(qJ(5));
t119 = cos(pkin(9));
t152 = t119 * t114;
t97 = -qJD(5) + t152;
t144 = t114 * t121 * t97;
t111 = qJDD(1) + qJDD(3);
t124 = cos(qJ(5));
t156 = t111 * t124;
t147 = qJD(5) * t114;
t93 = t121 * t117 * t147;
t130 = (t144 + t156) * t117 - t93;
t170 = t130 * t117;
t163 = t119 * pkin(4);
t164 = t117 * pkin(7);
t137 = -t163 - t164;
t166 = 2 * qJD(4);
t110 = t114 ^ 2;
t122 = sin(qJ(3));
t125 = cos(qJ(3));
t118 = sin(pkin(8));
t120 = cos(pkin(8));
t123 = sin(qJ(1));
t126 = cos(qJ(1));
t148 = t126 * g(2) + t123 * g(3);
t90 = qJDD(1) * pkin(1) + t148;
t127 = qJD(1) ^ 2;
t136 = -t123 * g(2) + t126 * g(3);
t91 = -t127 * pkin(1) - t136;
t140 = -t118 * t91 + t120 * t90;
t63 = qJDD(1) * pkin(2) + t140;
t161 = t118 * t90 + t120 * t91;
t64 = -t127 * pkin(2) + t161;
t48 = t122 * t63 + t125 * t64;
t46 = -t110 * pkin(3) + t111 * qJ(4) + t48;
t169 = t114 * t166 + t46;
t112 = t117 ^ 2;
t113 = t119 ^ 2;
t85 = (t112 + t113) * t110;
t149 = -g(1) + qJDD(2);
t105 = t119 * t149;
t30 = t169 * t117 - t105;
t31 = t117 * t149 + t169 * t119;
t14 = t117 * t30 + t119 * t31;
t157 = t111 * t121;
t168 = (t124 * t147 + t157) * t117;
t94 = t97 ^ 2;
t79 = t137 * t114;
t21 = -t105 + (t46 + (t166 + t79) * t114) * t117;
t22 = t79 * t152 + t31;
t47 = -t122 * t64 + t125 * t63;
t43 = -t111 * pkin(3) - t110 * qJ(4) + qJDD(4) - t47;
t33 = t137 * t111 + t43;
t11 = t121 * t22 - t124 * t33;
t12 = t121 * t33 + t124 * t22;
t7 = t121 * t11 + t124 * t12;
t3 = t117 * t21 + t119 * t7;
t6 = -t124 * t11 + t121 * t12;
t165 = -pkin(3) * t6 + qJ(4) * t3;
t162 = -pkin(3) * t43 + qJ(4) * t14;
t158 = t110 * t112;
t143 = t121 * t158;
t89 = t124 * t143;
t153 = t119 * t111;
t95 = -qJDD(5) + t153;
t68 = -t89 + t95;
t160 = t121 * t68;
t69 = -t89 - t95;
t159 = t124 * t69;
t155 = t114 * t124;
t154 = t117 * t111;
t151 = t122 * t111;
t150 = t125 * t111;
t77 = t119 * t85;
t145 = pkin(3) * t153 - qJ(4) * t77 - t119 * t43;
t142 = t124 ^ 2 * t158;
t67 = -t142 - t94;
t53 = -t121 * t67 + t124 * t68;
t32 = t119 * t53 + t170;
t52 = t124 * t67 + t160;
t141 = t117 * (-pkin(7) * t52 + t124 * t21) + t119 * (-pkin(4) * t52 + t12) - pkin(3) * t52 + qJ(4) * t32;
t96 = t121 ^ 2 * t158;
t70 = -t96 - t94;
t56 = -t121 * t69 + t124 * t70;
t75 = t97 * t117 * t155;
t58 = t75 - t168;
t37 = -t117 * t58 + t119 * t56;
t55 = t121 * t70 + t159;
t139 = t117 * (-pkin(7) * t55 + t121 * t21) + t119 * (-pkin(4) * t55 + t11) - pkin(3) * t55 + qJ(4) * t37;
t102 = t112 * t111;
t103 = t113 * t111;
t83 = t103 + t102;
t138 = pkin(3) * t85 + qJ(4) * t83 + t14;
t76 = t117 * t85;
t135 = -pkin(3) * t154 + qJ(4) * t76 + t117 * t43;
t86 = -t125 * t110 - t151;
t134 = t122 * t110 - t150;
t57 = t75 + t168;
t59 = -t93 + (-t144 + t156) * t117;
t45 = t121 * t59 - t124 * t57;
t71 = t96 + t142;
t25 = -t117 * t71 + t119 * t45;
t44 = -t121 * t57 - t124 * t59;
t133 = qJ(4) * t25 - t117 * t6 + (-pkin(3) + t137) * t44;
t92 = 0.2e1 * t117 * t153;
t88 = t119 * t95;
t72 = -t96 + t142;
t66 = t119 * t150 - t122 * t77;
t65 = -t117 * t150 + t122 * t76;
t62 = t122 * t83 + t125 * t85;
t50 = (-t119 * t143 + t170) * t124;
t49 = t119 * t89 + (t157 + (qJD(5) - t97) * t155) * t121 * t112;
t36 = t117 * (t124 * (t96 - t94) + t160) + t119 * t57;
t35 = t117 * (t159 - t121 * (t94 - t142)) - t119 * t59;
t24 = t117 * (-t121 * t130 + t124 * t58) - t119 * t72;
t20 = t122 * t37 - t125 * t55;
t19 = t122 * t48 + t125 * t47;
t18 = t122 * t32 - t125 * t52;
t15 = t122 * t25 - t125 * t44;
t8 = t122 * t14 - t125 * t43;
t1 = t122 * t3 - t125 * t6;
t2 = [0, 0, 0, 0, 0, qJDD(1), t148, t136, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t120 * qJDD(1) - t118 * t127) + t140, pkin(1) * (-t118 * qJDD(1) - t120 * t127) - t161, 0, pkin(1) * (t118 * t161 + t120 * t140), 0, 0, 0, 0, 0, t111, pkin(1) * (t118 * t86 - t120 * t134) - pkin(2) * t134 + t47, pkin(1) * (t118 * t134 + t120 * t86) + pkin(2) * t86 - t48, 0, pkin(1) * (t118 * (-t122 * t47 + t125 * t48) + t120 * t19) + pkin(2) * t19, t102, t92, 0, t103, 0, 0, pkin(1) * (t118 * (-t119 * t151 - t125 * t77) + t120 * t66) + pkin(2) * t66 + t145, pkin(1) * (t118 * (t117 * t151 + t125 * t76) + t120 * t65) + pkin(2) * t65 + t135, pkin(1) * (t118 * (-t122 * t85 + t125 * t83) + t120 * t62) + pkin(2) * t62 + t138, pkin(1) * (t118 * (t122 * t43 + t125 * t14) + t120 * t8) + pkin(2) * t8 + t162, t50, t24, t35, t49, t36, t88, pkin(1) * (t118 * (t122 * t55 + t125 * t37) + t120 * t20) + pkin(2) * t20 + t139, pkin(1) * (t118 * (t122 * t52 + t125 * t32) + t120 * t18) + pkin(2) * t18 + t141, pkin(1) * (t118 * (t122 * t44 + t125 * t25) + t120 * t15) + pkin(2) * t15 + t133, pkin(1) * (t118 * (t122 * t6 + t125 * t3) + t120 * t1) + pkin(2) * t1 - t6 * t164 - t6 * t163 + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117 * t31 - t119 * t30, 0, 0, 0, 0, 0, 0, t117 * t56 + t119 * t58, t117 * t53 - t119 * t130, t117 * t45 + t119 * t71, t117 * t7 - t119 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t47, -t48, 0, 0, t102, t92, 0, t103, 0, 0, t145, t135, t138, t162, t50, t24, t35, t49, t36, t88, t139, t141, t133, t137 * t6 + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t154, -t85, t43, 0, 0, 0, 0, 0, 0, t55, t52, t44, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t72, t59, -t89, -t57, -t95, -t11, -t12, 0, 0;];
tauJ_reg = t2;
