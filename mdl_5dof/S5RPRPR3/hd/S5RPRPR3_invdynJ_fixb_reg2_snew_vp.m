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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:36:22
% EndTime: 2020-01-03 11:36:29
% DurationCPUTime: 1.05s
% Computational Cost: add. (4082->165), mult. (6022->245), div. (0->0), fcn. (3535->10), ass. (0->123)
t117 = cos(pkin(9));
t162 = t117 * pkin(4);
t115 = sin(pkin(9));
t163 = t115 * pkin(7);
t137 = -t162 - t163;
t112 = (qJD(1) + qJD(3));
t165 = 2 * qJD(4);
t108 = t112 ^ 2;
t109 = qJDD(1) + qJDD(3);
t120 = sin(qJ(3));
t123 = cos(qJ(3));
t116 = sin(pkin(8));
t118 = cos(pkin(8));
t121 = sin(qJ(1));
t124 = cos(qJ(1));
t136 = -t124 * g(2) - t121 * g(3);
t90 = qJDD(1) * pkin(1) + t136;
t125 = qJD(1) ^ 2;
t135 = t121 * g(2) - t124 * g(3);
t91 = -t125 * pkin(1) - t135;
t140 = -t116 * t91 + t118 * t90;
t63 = qJDD(1) * pkin(2) + t140;
t160 = t116 * t90 + t118 * t91;
t64 = -t125 * pkin(2) + t160;
t48 = t120 * t63 + t123 * t64;
t46 = -t108 * pkin(3) + t109 * qJ(4) + t48;
t168 = (t112 * t165) + t46;
t110 = t115 ^ 2;
t111 = t117 ^ 2;
t85 = (t110 + t111) * t108;
t147 = -g(1) + qJDD(2);
t105 = t117 * t147;
t30 = t168 * t115 - t105;
t31 = t115 * t147 + t168 * t117;
t14 = t115 * t30 + t117 * t31;
t122 = cos(qJ(5));
t146 = qJD(5) * t112;
t119 = sin(qJ(5));
t155 = t109 * t119;
t167 = (t122 * t146 + t155) * t115;
t150 = t117 * t112;
t97 = -qJD(5) + t150;
t159 = t112 * t97;
t143 = t119 * t159;
t154 = t109 * t122;
t93 = t119 * t115 * t146;
t128 = (t143 + t154) * t115 - t93;
t94 = t97 ^ 2;
t79 = t137 * t112;
t21 = -t105 + (t46 + (t165 + t79) * t112) * t115;
t22 = t79 * t150 + t31;
t47 = -t120 * t64 + t123 * t63;
t43 = -t109 * pkin(3) - t108 * qJ(4) + qJDD(4) - t47;
t33 = t137 * t109 + t43;
t11 = t119 * t22 - t122 * t33;
t12 = t119 * t33 + t122 * t22;
t7 = t119 * t11 + t122 * t12;
t3 = t115 * t21 + t117 * t7;
t6 = -t122 * t11 + t119 * t12;
t164 = -pkin(3) * t6 + qJ(4) * t3;
t161 = -pkin(3) * t43 + qJ(4) * t14;
t156 = t108 * t110;
t89 = t119 * t122 * t156;
t151 = t117 * t109;
t95 = -qJDD(5) + t151;
t68 = -t89 + t95;
t158 = t119 * t68;
t69 = -t89 - t95;
t157 = t122 * t69;
t153 = t115 * t109;
t152 = t115 * t122;
t149 = t120 * t109;
t148 = t123 * t109;
t77 = t117 * t85;
t144 = pkin(3) * t151 - qJ(4) * t77 - t117 * t43;
t142 = t122 ^ 2 * t156;
t67 = -t142 - t94;
t53 = -t119 * t67 + t122 * t68;
t32 = t115 * t128 + t117 * t53;
t52 = t122 * t67 + t158;
t141 = t115 * (-pkin(7) * t52 + t122 * t21) + t117 * (-pkin(4) * t52 + t12) - pkin(3) * t52 + qJ(4) * t32;
t96 = t119 ^ 2 * t156;
t70 = -t96 - t94;
t56 = -t119 * t69 + t122 * t70;
t75 = t152 * t159;
t58 = t75 - t167;
t37 = -t115 * t58 + t117 * t56;
t55 = t119 * t70 + t157;
t139 = t115 * (-pkin(7) * t55 + t119 * t21) + t117 * (-pkin(4) * t55 + t11) - pkin(3) * t55 + qJ(4) * t37;
t102 = t110 * t109;
t103 = t111 * t109;
t83 = t103 + t102;
t138 = pkin(3) * t85 + qJ(4) * t83 + t14;
t134 = t117 * t89;
t76 = t115 * t85;
t133 = -pkin(3) * t153 + qJ(4) * t76 + t115 * t43;
t86 = -t123 * t108 - t149;
t132 = t120 * t108 - t148;
t57 = t75 + t167;
t59 = -t93 + (-t143 + t154) * t115;
t45 = t119 * t59 - t122 * t57;
t71 = t96 + t142;
t25 = -t115 * t71 + t117 * t45;
t44 = -t119 * t57 - t122 * t59;
t131 = qJ(4) * t25 - t115 * t6 + (-pkin(3) + t137) * t44;
t92 = 0.2e1 * t115 * t151;
t88 = t117 * t95;
t72 = -t96 + t142;
t66 = t117 * t148 - t120 * t77;
t65 = -t115 * t148 + t120 * t76;
t62 = t120 * t83 + t123 * t85;
t50 = t128 * t152 - t134;
t49 = t134 + (t155 + (qJD(5) - t97) * t122 * t112) * t119 * t110;
t36 = -t117 * t59 + t115 * (t157 - t119 * (t94 - t142));
t35 = t117 * t57 + t115 * (t122 * (t96 - t94) + t158);
t24 = -t117 * t72 + t115 * (-t119 * t128 + t122 * t58);
t20 = t120 * t37 - t123 * t55;
t19 = t120 * t48 + t123 * t47;
t18 = t120 * t32 - t123 * t52;
t15 = t120 * t25 - t123 * t44;
t8 = t120 * t14 - t123 * t43;
t1 = t120 * t3 - t123 * t6;
t2 = [0, 0, 0, 0, 0, qJDD(1), t136, t135, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t118 * qJDD(1) - t116 * t125) + t140, pkin(1) * (-t116 * qJDD(1) - t118 * t125) - t160, 0, pkin(1) * (t116 * t160 + t118 * t140), 0, 0, 0, 0, 0, t109, pkin(1) * (t116 * t86 - t118 * t132) - pkin(2) * t132 + t47, pkin(1) * (t116 * t132 + t118 * t86) + pkin(2) * t86 - t48, 0, pkin(1) * (t116 * (-t120 * t47 + t123 * t48) + t118 * t19) + pkin(2) * t19, t102, t92, 0, t103, 0, 0, pkin(1) * (t116 * (-t117 * t149 - t123 * t77) + t118 * t66) + pkin(2) * t66 + t144, pkin(1) * (t116 * (t115 * t149 + t123 * t76) + t118 * t65) + pkin(2) * t65 + t133, pkin(1) * (t116 * (-t120 * t85 + t123 * t83) + t118 * t62) + pkin(2) * t62 + t138, pkin(1) * (t116 * (t120 * t43 + t123 * t14) + t118 * t8) + pkin(2) * t8 + t161, t50, t24, t36, t49, t35, t88, pkin(1) * (t116 * (t120 * t55 + t123 * t37) + t118 * t20) + pkin(2) * t20 + t139, pkin(1) * (t116 * (t120 * t52 + t123 * t32) + t118 * t18) + pkin(2) * t18 + t141, pkin(1) * (t116 * (t120 * t44 + t123 * t25) + t118 * t15) + pkin(2) * t15 + t131, pkin(1) * (t116 * (t120 * t6 + t123 * t3) + t118 * t1) + pkin(2) * t1 - t6 * t162 - t6 * t163 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115 * t31 - t117 * t30, 0, 0, 0, 0, 0, 0, t115 * t56 + t117 * t58, t115 * t53 - t117 * t128, t115 * t45 + t117 * t71, t115 * t7 - t117 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, t47, -t48, 0, 0, t102, t92, 0, t103, 0, 0, t144, t133, t138, t161, t50, t24, t36, t49, t35, t88, t139, t141, t131, t137 * t6 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t153, -t85, t43, 0, 0, 0, 0, 0, 0, t55, t52, t44, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t72, t59, -t89, -t57, -t95, -t11, -t12, 0, 0;];
tauJ_reg = t2;
