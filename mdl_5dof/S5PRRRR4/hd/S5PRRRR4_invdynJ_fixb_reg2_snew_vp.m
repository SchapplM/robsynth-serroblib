% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:59
% EndTime: 2019-12-05 17:08:03
% DurationCPUTime: 1.00s
% Computational Cost: add. (4370->171), mult. (5962->235), div. (0->0), fcn. (4085->10), ass. (0->126)
t129 = sin(qJ(5));
t120 = qJDD(4) + qJDD(5);
t123 = qJD(2) + qJD(3);
t133 = cos(qJ(5));
t134 = cos(qJ(4));
t130 = sin(qJ(4));
t159 = t123 * t130;
t86 = -t133 * t134 * t123 + t129 * t159;
t88 = (t134 * t129 + t130 * t133) * t123;
t67 = t88 * t86;
t170 = -t67 + t120;
t173 = t129 * t170;
t172 = t133 * t170;
t126 = -g(3) + qJDD(1);
t119 = t123 ^ 2;
t121 = qJDD(2) + qJDD(3);
t131 = sin(qJ(3));
t135 = cos(qJ(3));
t127 = sin(pkin(9));
t128 = cos(pkin(9));
t104 = t127 * g(1) - t128 * g(2);
t105 = -t128 * g(1) - t127 * g(2);
t132 = sin(qJ(2));
t136 = cos(qJ(2));
t143 = t136 * t104 - t132 * t105;
t139 = qJDD(2) * pkin(2) + t143;
t142 = -t132 * t104 - t136 * t105;
t74 = -qJD(2) ^ 2 * pkin(2) - t142;
t57 = t131 * t139 + t135 * t74;
t55 = -t119 * pkin(3) + t121 * pkin(7) + t57;
t39 = -t134 * t126 + t130 * t55;
t40 = t130 * t126 + t134 * t55;
t23 = t130 * t39 + t134 * t40;
t155 = qJD(4) * t123;
t149 = t134 * t155;
t157 = t130 * t121;
t95 = t149 + t157;
t112 = t134 * t121;
t150 = t130 * t155;
t96 = t112 - t150;
t60 = -t86 * qJD(5) + t129 * t96 + t133 * t95;
t122 = qJD(4) + qJD(5);
t83 = t122 * t86;
t171 = t60 - t83;
t169 = -t39 + (-t95 + t149) * pkin(8);
t84 = t86 ^ 2;
t85 = t88 ^ 2;
t118 = t122 ^ 2;
t108 = t130 * t119 * t134;
t154 = qJDD(4) + t108;
t138 = t154 * pkin(4) + t169;
t103 = qJD(4) * pkin(4) - pkin(8) * t159;
t125 = t134 ^ 2;
t114 = t125 * t119;
t31 = -pkin(4) * t114 + t96 * pkin(8) - qJD(4) * t103 + t40;
t14 = t129 * t31 - t133 * t138;
t164 = t133 * t31;
t15 = t129 * t138 + t164;
t6 = t129 * t15 - t133 * t14;
t168 = t130 * t6;
t56 = -t131 * t74 + t135 * t139;
t54 = -t121 * pkin(3) - t119 * pkin(7) - t56;
t167 = -pkin(3) * t54 + pkin(7) * t23;
t33 = -t96 * pkin(4) - pkin(8) * t114 + t103 * t159 + t54;
t166 = t129 * t33;
t64 = t67 + t120;
t165 = t129 * t64;
t163 = t133 * t33;
t162 = t133 * t64;
t161 = t122 * t129;
t160 = t122 * t133;
t158 = t130 * t154;
t102 = qJDD(4) - t108;
t156 = t134 * t102;
t124 = t130 ^ 2;
t113 = t124 * t119;
t137 = qJD(4) ^ 2;
t106 = -t113 - t137;
t78 = -t130 * t106 - t156;
t94 = 0.2e1 * t149 + t157;
t153 = -pkin(3) * t94 + pkin(7) * t78 + t130 * t54;
t107 = -t114 - t137;
t77 = t134 * t107 - t158;
t97 = t112 - 0.2e1 * t150;
t152 = pkin(3) * t97 + pkin(7) * t77 - t134 * t54;
t148 = t129 * t95 - t133 * t96;
t140 = (-qJD(5) + t122) * t88 - t148;
t48 = t60 + t83;
t24 = t129 * t140 - t133 * t48;
t25 = t129 * t48 + t133 * t140;
t10 = -t130 * t24 + t134 * t25;
t61 = -t84 - t85;
t7 = t129 * t14 + t133 * t15;
t151 = t130 * (-pkin(8) * t24 - t6) + t134 * (-pkin(4) * t61 + pkin(8) * t25 + t7) - pkin(3) * t61 + pkin(7) * t10;
t62 = -t118 - t84;
t35 = t129 * t62 + t172;
t36 = t133 * t62 - t173;
t19 = -t130 * t35 + t134 * t36;
t43 = (qJD(5) + t122) * t88 + t148;
t147 = t130 * (-pkin(8) * t35 + t166) + t134 * (-pkin(4) * t43 + pkin(8) * t36 - t163) - pkin(3) * t43 + pkin(7) * t19;
t79 = -t85 - t118;
t51 = t133 * t79 - t165;
t52 = -t129 * t79 - t162;
t27 = -t130 * t51 + t134 * t52;
t146 = t130 * (-pkin(8) * t51 + t163) + t134 * (-pkin(4) * t171 + pkin(8) * t52 + t166) - pkin(3) * t171 + pkin(7) * t27;
t100 = t113 + t114;
t99 = (t124 + t125) * t121;
t145 = pkin(3) * t100 + pkin(7) * t99 + t23;
t2 = t134 * t7 - t168;
t141 = pkin(7) * t2 - pkin(8) * t168 - pkin(3) * t33 + t134 * (-pkin(4) * t33 + pkin(8) * t7);
t81 = -t85 + t118;
t80 = t84 - t118;
t76 = t158 + t134 * (-t113 + t137);
t75 = t130 * (t114 - t137) + t156;
t71 = (t95 + t149) * t130;
t70 = (t96 - t150) * t134;
t68 = t130 * t97 + t134 * t94;
t66 = t85 - t84;
t59 = -t88 * qJD(5) - t148;
t34 = (t130 * (t129 * t88 - t133 * t86) + t134 * (-t129 * t86 - t133 * t88)) * t122;
t29 = t130 * (t133 * t80 - t165) + t134 * (t129 * t80 + t162);
t28 = t130 * (-t129 * t81 + t172) + t134 * (t133 * t81 + t173);
t21 = t130 * (t133 * t60 - t88 * t161) + t134 * (t129 * t60 + t88 * t160);
t20 = t130 * (-t129 * t59 + t86 * t160) + t134 * (t133 * t59 + t86 * t161);
t9 = t130 * (-t129 * t171 - t133 * t43) + t134 * (-t129 * t43 + t133 * t171);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0, 0, 0, 0, 0, 0, t130 * t107 + t134 * t154, -t130 * t102 + t134 * t106, 0, t130 * t40 - t134 * t39, 0, 0, 0, 0, 0, 0, t130 * t36 + t134 * t35, t130 * t52 + t134 * t51, t130 * t25 + t134 * t24, t130 * t7 + t134 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t143, t142, 0, 0, 0, 0, 0, 0, 0, t121, pkin(2) * (-t131 * t119 + t135 * t121) + t56, pkin(2) * (-t135 * t119 - t131 * t121) - t57, 0, pkin(2) * (t131 * t57 + t135 * t56), t71, t68, t76, t70, t75, 0, pkin(2) * (t131 * t77 + t135 * t97) + t152, pkin(2) * (t131 * t78 - t135 * t94) + t153, pkin(2) * (t135 * t100 + t131 * t99) + t145, pkin(2) * (t131 * t23 - t135 * t54) + t167, t21, t9, t28, t20, t29, t34, pkin(2) * (t131 * t19 - t135 * t43) + t147, pkin(2) * (t131 * t27 - t135 * t171) + t146, pkin(2) * (t131 * t10 - t135 * t61) + t151, pkin(2) * (t131 * t2 - t135 * t33) + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t56, -t57, 0, 0, t71, t68, t76, t70, t75, 0, t152, t153, t145, t167, t21, t9, t28, t20, t29, t34, t147, t146, t151, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, t113 - t114, t157, t108, t112, qJDD(4), -t39, -t40, 0, 0, t67, t66, t48, -t67, t140, t120, pkin(4) * t35 - t14, -t164 - t129 * t169 + (-t129 * t154 + t51) * pkin(4), pkin(4) * t24, pkin(4) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t66, t48, -t67, t140, t120, -t14, -t15, 0, 0;];
tauJ_reg = t1;
