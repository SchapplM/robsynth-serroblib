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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:20:55
% EndTime: 2022-01-23 09:20:58
% DurationCPUTime: 1.02s
% Computational Cost: add. (4082->165), mult. (6022->245), div. (0->0), fcn. (3535->10), ass. (0->123)
t116 = sin(pkin(9));
t113 = qJD(1) + qJD(3);
t120 = sin(qJ(5));
t118 = cos(pkin(9));
t151 = t118 * t113;
t97 = -qJD(5) + t151;
t144 = t113 * t120 * t97;
t110 = qJDD(1) + qJDD(3);
t123 = cos(qJ(5));
t156 = t110 * t123;
t147 = qJD(5) * t113;
t93 = t120 * t116 * t147;
t129 = (t144 + t156) * t116 - t93;
t169 = t129 * t116;
t162 = t118 * pkin(4);
t163 = t116 * pkin(7);
t136 = -t162 - t163;
t165 = 2 * qJD(4);
t109 = t113 ^ 2;
t121 = sin(qJ(3));
t124 = cos(qJ(3));
t117 = sin(pkin(8));
t119 = cos(pkin(8));
t122 = sin(qJ(1));
t125 = cos(qJ(1));
t140 = t122 * g(1) - t125 * g(2);
t90 = qJDD(1) * pkin(1) + t140;
t126 = qJD(1) ^ 2;
t135 = t125 * g(1) + t122 * g(2);
t91 = -t126 * pkin(1) - t135;
t139 = -t117 * t91 + t119 * t90;
t63 = qJDD(1) * pkin(2) + t139;
t160 = t117 * t90 + t119 * t91;
t64 = -t126 * pkin(2) + t160;
t48 = t121 * t63 + t124 * t64;
t46 = -t109 * pkin(3) + t110 * qJ(4) + t48;
t168 = t113 * t165 + t46;
t111 = t116 ^ 2;
t112 = t118 ^ 2;
t85 = (t111 + t112) * t109;
t148 = -g(3) + qJDD(2);
t105 = t118 * t148;
t30 = t168 * t116 - t105;
t31 = t116 * t148 + t168 * t118;
t14 = t116 * t30 + t118 * t31;
t157 = t110 * t120;
t167 = (t123 * t147 + t157) * t116;
t94 = t97 ^ 2;
t79 = t136 * t113;
t21 = -t105 + (t46 + (t165 + t79) * t113) * t116;
t22 = t79 * t151 + t31;
t47 = -t121 * t64 + t124 * t63;
t43 = -t110 * pkin(3) - t109 * qJ(4) + qJDD(4) - t47;
t33 = t136 * t110 + t43;
t11 = t120 * t22 - t123 * t33;
t12 = t120 * t33 + t123 * t22;
t7 = t120 * t11 + t123 * t12;
t3 = t116 * t21 + t118 * t7;
t6 = -t123 * t11 + t120 * t12;
t164 = -pkin(3) * t6 + qJ(4) * t3;
t161 = -pkin(3) * t43 + qJ(4) * t14;
t155 = t111 * t109;
t143 = t120 * t155;
t89 = t123 * t143;
t152 = t118 * t110;
t95 = -qJDD(5) + t152;
t68 = -t89 + t95;
t159 = t120 * t68;
t69 = -t89 - t95;
t158 = t123 * t69;
t154 = t113 * t123;
t153 = t116 * t110;
t150 = t121 * t110;
t149 = t124 * t110;
t77 = t118 * t85;
t145 = pkin(3) * t152 - qJ(4) * t77 - t118 * t43;
t142 = t123 ^ 2 * t155;
t67 = -t142 - t94;
t53 = -t120 * t67 + t123 * t68;
t32 = t118 * t53 + t169;
t52 = t123 * t67 + t159;
t141 = t116 * (-pkin(7) * t52 + t123 * t21) + t118 * (-pkin(4) * t52 + t12) - pkin(3) * t52 + qJ(4) * t32;
t96 = t120 ^ 2 * t155;
t70 = -t96 - t94;
t56 = -t120 * t69 + t123 * t70;
t75 = t97 * t116 * t154;
t58 = t75 - t167;
t37 = -t116 * t58 + t118 * t56;
t55 = t120 * t70 + t158;
t138 = t116 * (-pkin(7) * t55 + t120 * t21) + t118 * (-pkin(4) * t55 + t11) - pkin(3) * t55 + qJ(4) * t37;
t102 = t111 * t110;
t103 = t112 * t110;
t83 = t103 + t102;
t137 = pkin(3) * t85 + qJ(4) * t83 + t14;
t76 = t116 * t85;
t134 = -pkin(3) * t153 + qJ(4) * t76 + t116 * t43;
t86 = -t124 * t109 - t150;
t133 = t121 * t109 - t149;
t57 = t75 + t167;
t59 = -t93 + (-t144 + t156) * t116;
t45 = t120 * t59 - t123 * t57;
t71 = t96 + t142;
t25 = -t116 * t71 + t118 * t45;
t44 = -t120 * t57 - t123 * t59;
t132 = qJ(4) * t25 - t116 * t6 + (-pkin(3) + t136) * t44;
t92 = 0.2e1 * t116 * t152;
t88 = t118 * t95;
t72 = -t96 + t142;
t66 = t118 * t149 - t121 * t77;
t65 = -t116 * t149 + t121 * t76;
t62 = t121 * t83 + t124 * t85;
t50 = (-t118 * t143 + t169) * t123;
t49 = t118 * t89 + (t157 + (qJD(5) - t97) * t154) * t120 * t111;
t36 = t116 * (t123 * (t96 - t94) + t159) + t118 * t57;
t35 = t116 * (t158 - t120 * (t94 - t142)) - t118 * t59;
t24 = t116 * (-t120 * t129 + t123 * t58) - t118 * t72;
t20 = t121 * t37 - t124 * t55;
t19 = t121 * t48 + t124 * t47;
t18 = t121 * t32 - t124 * t52;
t15 = t121 * t25 - t124 * t44;
t8 = t121 * t14 - t124 * t43;
t1 = t121 * t3 - t124 * t6;
t2 = [0, 0, 0, 0, 0, qJDD(1), t140, t135, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t119 * qJDD(1) - t117 * t126) + t139, pkin(1) * (-t117 * qJDD(1) - t119 * t126) - t160, 0, pkin(1) * (t117 * t160 + t119 * t139), 0, 0, 0, 0, 0, t110, pkin(1) * (t117 * t86 - t119 * t133) - pkin(2) * t133 + t47, pkin(1) * (t117 * t133 + t119 * t86) + pkin(2) * t86 - t48, 0, pkin(1) * (t117 * (-t121 * t47 + t124 * t48) + t119 * t19) + pkin(2) * t19, t102, t92, 0, t103, 0, 0, pkin(1) * (t117 * (-t118 * t150 - t124 * t77) + t119 * t66) + pkin(2) * t66 + t145, pkin(1) * (t117 * (t116 * t150 + t124 * t76) + t119 * t65) + pkin(2) * t65 + t134, pkin(1) * (t117 * (-t121 * t85 + t124 * t83) + t119 * t62) + pkin(2) * t62 + t137, pkin(1) * (t117 * (t121 * t43 + t124 * t14) + t119 * t8) + pkin(2) * t8 + t161, t50, t24, t35, t49, t36, t88, pkin(1) * (t117 * (t121 * t55 + t124 * t37) + t119 * t20) + pkin(2) * t20 + t138, pkin(1) * (t117 * (t121 * t52 + t124 * t32) + t119 * t18) + pkin(2) * t18 + t141, pkin(1) * (t117 * (t121 * t44 + t124 * t25) + t119 * t15) + pkin(2) * t15 + t132, pkin(1) * (t117 * (t121 * t6 + t124 * t3) + t119 * t1) + pkin(2) * t1 - t6 * t163 - t6 * t162 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116 * t31 - t118 * t30, 0, 0, 0, 0, 0, 0, t116 * t56 + t118 * t58, t116 * t53 - t118 * t129, t116 * t45 + t118 * t71, t116 * t7 - t118 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t47, -t48, 0, 0, t102, t92, 0, t103, 0, 0, t145, t134, t137, t161, t50, t24, t35, t49, t36, t88, t138, t141, t132, t136 * t6 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, t153, -t85, t43, 0, 0, 0, 0, 0, 0, t55, t52, t44, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t72, t59, -t89, -t57, -t95, -t11, -t12, 0, 0;];
tauJ_reg = t2;
