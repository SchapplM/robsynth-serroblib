% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:46
% EndTime: 2019-12-05 15:19:51
% DurationCPUTime: 1.45s
% Computational Cost: add. (4808->210), mult. (8808->325), div. (0->0), fcn. (7272->14), ass. (0->150)
t106 = sin(pkin(6));
t110 = cos(pkin(6));
t104 = sin(pkin(11));
t108 = cos(pkin(11));
t101 = -g(3) + qJDD(1);
t107 = sin(pkin(5));
t111 = cos(pkin(5));
t105 = sin(pkin(10));
t109 = cos(pkin(10));
t83 = g(1) * t105 - g(2) * t109;
t132 = t101 * t107 + t111 * t83;
t144 = -g(1) * t109 - g(2) * t105;
t124 = -t104 * t144 + t108 * t132;
t126 = t101 * t111 - t107 * t83 + qJDD(2);
t173 = t106 * t126 + t110 * t124;
t114 = sin(qJ(5));
t117 = cos(qJ(5));
t115 = sin(qJ(4));
t153 = qJD(3) * t115;
t73 = -qJD(4) * t117 + t114 * t153;
t75 = qJD(4) * t114 + t117 * t153;
t169 = t75 * t73;
t151 = qJD(3) * qJD(4);
t93 = t115 * t151;
t118 = cos(qJ(4));
t95 = t118 * qJDD(3);
t79 = t95 - t93;
t72 = -qJDD(5) + t79;
t129 = -t72 - t169;
t172 = t114 * t129;
t171 = t117 * t129;
t116 = sin(qJ(3));
t119 = cos(qJ(3));
t49 = t104 * t132 + t108 * t144;
t146 = t116 * t49 - t119 * t173;
t32 = t116 * t173 + t119 * t49;
t170 = t116 * t32 - t119 * t146;
t148 = t118 * t151;
t150 = t115 * qJDD(3);
t78 = t148 + t150;
t147 = -qJDD(4) * t117 + t114 * t78;
t152 = t118 * qJD(3);
t90 = -qJD(5) + t152;
t42 = (qJD(5) + t90) * t75 + t147;
t70 = t73 ^ 2;
t71 = t75 ^ 2;
t88 = t90 ^ 2;
t120 = qJD(4) ^ 2;
t121 = qJD(3) ^ 2;
t30 = -t121 * pkin(3) + qJDD(3) * pkin(8) + t32;
t122 = -t106 * t124 + t110 * t126;
t41 = t118 * t122;
t145 = -pkin(4) * t118 - pkin(9) * t115;
t76 = t145 * qJD(3);
t16 = -qJDD(4) * pkin(4) - t120 * pkin(9) - t41 + (qJD(3) * t76 + t30) * t115;
t168 = t114 * t16;
t51 = t72 - t169;
t167 = t114 * t51;
t166 = t114 * t90;
t89 = t115 * t121 * t118;
t84 = qJDD(4) + t89;
t165 = t115 * t84;
t163 = t117 * t16;
t162 = t117 * t51;
t161 = t117 * t90;
t85 = qJDD(4) - t89;
t160 = t118 * t85;
t158 = t107 * t104;
t157 = t107 * t108;
t156 = t108 * t110;
t155 = qJD(5) - t90;
t149 = t118 * t169;
t20 = t115 * t122 + t118 * t30;
t17 = -t120 * pkin(4) + qJDD(4) * pkin(9) + t152 * t76 + t20;
t141 = t78 + t148;
t142 = -t79 + t93;
t29 = -qJDD(3) * pkin(3) - pkin(8) * t121 + t146;
t22 = pkin(4) * t142 - pkin(9) * t141 + t29;
t11 = t114 * t17 - t117 * t22;
t12 = t114 * t22 + t117 * t17;
t6 = t11 * t114 + t117 * t12;
t19 = t115 * t30 - t41;
t9 = t115 * t19 + t118 * t20;
t139 = t11 * t117 - t114 * t12;
t3 = t115 * t16 + t118 * t6;
t143 = t116 * t3 + t119 * t139;
t140 = t116 * t9 - t119 * t29;
t131 = -t114 * qJDD(4) - t117 * t78;
t55 = -qJD(5) * t73 - t131;
t67 = t73 * t90;
t46 = t55 - t67;
t34 = t114 * t46 - t117 * t42;
t50 = t70 + t71;
t24 = -t115 * t50 + t118 * t34;
t33 = -t114 * t42 - t117 * t46;
t138 = t116 * t24 - t119 * t33;
t56 = -t88 - t70;
t36 = t117 * t56 - t172;
t43 = -t155 * t75 - t147;
t26 = -t115 * t43 + t118 * t36;
t35 = t114 * t56 + t171;
t137 = t116 * t26 - t119 * t35;
t59 = -t71 - t88;
t40 = -t114 * t59 + t162;
t47 = t155 * t73 + t131;
t28 = -t115 * t47 + t118 * t40;
t39 = t117 * t59 + t167;
t136 = t116 * t28 - t119 * t39;
t100 = t118 ^ 2;
t98 = t100 * t121;
t87 = -t98 - t120;
t63 = t118 * t87 - t165;
t80 = t95 - 0.2e1 * t93;
t135 = t116 * t63 + t119 * t80;
t99 = t115 ^ 2;
t96 = t99 * t121;
t86 = -t96 - t120;
t64 = -t115 * t86 - t160;
t77 = 0.2e1 * t148 + t150;
t134 = t116 * t64 - t119 * t77;
t81 = (t100 + t99) * qJDD(3);
t82 = t96 + t98;
t133 = t116 * t81 + t119 * t82;
t130 = -pkin(3) + t145;
t128 = qJDD(3) * t119 - t116 * t121;
t127 = -qJDD(3) * t116 - t119 * t121;
t69 = t128 * t106;
t68 = t127 * t106;
t66 = -t71 + t88;
t65 = t70 - t88;
t62 = -t115 * t85 + t118 * t86;
t61 = t115 * t87 + t118 * t84;
t60 = t71 - t70;
t57 = t133 * t106;
t54 = -qJD(5) * t75 - t147;
t45 = t55 + t67;
t38 = t106 * t134 + t110 * t62;
t37 = t106 * t135 + t110 * t61;
t27 = t115 * t40 + t118 * t47;
t25 = t115 * t36 + t118 * t43;
t23 = t115 * t34 + t118 * t50;
t15 = t106 * t170 + t110 * t122;
t14 = t106 * t136 + t110 * t27;
t13 = t106 * t137 + t110 * t25;
t8 = t115 * t20 - t118 * t19;
t7 = t106 * t138 + t110 * t23;
t4 = t106 * t140 + t110 * t8;
t2 = t115 * t6 - t118 * t16;
t1 = t106 * t143 + t110 * t2;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111 * t126 + t124 * t157 + t158 * t49, 0, 0, 0, 0, 0, 0, t111 * t69 + (t104 * t127 + t128 * t156) * t107, t111 * t68 + (-t104 * t128 + t127 * t156) * t107, 0, (t116 * t146 + t119 * t32) * t158 + (-t106 * t122 + t110 * t170) * t157 + t111 * t15, 0, 0, 0, 0, 0, 0, t111 * t37 + (t104 * (-t116 * t80 + t119 * t63) + t108 * (-t106 * t61 + t110 * t135)) * t107, t111 * t38 + (t104 * (t116 * t77 + t119 * t64) + t108 * (-t106 * t62 + t110 * t134)) * t107, t111 * t57 + (t104 * (-t116 * t82 + t119 * t81) + t133 * t156) * t107, t111 * t4 + (t104 * (t116 * t29 + t119 * t9) + t108 * (-t106 * t8 + t110 * t140)) * t107, 0, 0, 0, 0, 0, 0, t111 * t13 + (t104 * (t116 * t35 + t119 * t26) + t108 * (-t106 * t25 + t110 * t137)) * t107, t111 * t14 + (t104 * (t116 * t39 + t119 * t28) + t108 * (-t106 * t27 + t110 * t136)) * t107, t111 * t7 + (t104 * (t116 * t33 + t119 * t24) + t108 * (-t106 * t23 + t110 * t138)) * t107, t111 * t1 + (t104 * (-t116 * t139 + t119 * t3) + t108 * (-t106 * t2 + t110 * t143)) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0, 0, 0, 0, 0, 0, t69, t68, 0, t15, 0, 0, 0, 0, 0, 0, t37, t38, t57, t4, 0, 0, 0, 0, 0, 0, t13, t14, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t146, -t32, 0, 0, t141 * t115, t115 * t80 + t118 * t77, t165 + t118 * (-t96 + t120), -t142 * t118, t115 * (t98 - t120) + t160, 0, pkin(3) * t80 + pkin(8) * t63 - t118 * t29, -pkin(3) * t77 + pkin(8) * t64 + t115 * t29, pkin(3) * t82 + pkin(8) * t81 + t9, -pkin(3) * t29 + pkin(8) * t9, t115 * (t117 * t55 + t166 * t75) - t149, t115 * (-t114 * t45 + t117 * t43) - t118 * t60, t115 * (-t114 * t66 + t171) - t118 * t46, t115 * (-t114 * t54 - t161 * t73) + t149, t115 * (t117 * t65 + t167) + t118 * t42, t118 * t72 + t115 * (-t114 * t75 + t117 * t73) * t90, t115 * (-pkin(9) * t35 + t168) + t118 * (-pkin(4) * t35 + t11) - pkin(3) * t35 + pkin(8) * t26, t115 * (-pkin(9) * t39 + t163) + t118 * (-pkin(4) * t39 + t12) - pkin(3) * t39 + pkin(8) * t28, pkin(8) * t24 + t115 * t139 + t130 * t33, pkin(8) * t3 - t130 * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t96 - t98, t150, t89, t95, qJDD(4), -t19, -t20, 0, 0, t114 * t55 - t161 * t75, t114 * t43 + t117 * t45, t117 * t66 + t172, t117 * t54 - t166 * t73, t114 * t65 - t162, (t114 * t73 + t117 * t75) * t90, pkin(4) * t43 + pkin(9) * t36 - t163, pkin(4) * t47 + pkin(9) * t40 + t168, pkin(4) * t50 + pkin(9) * t34 + t6, -pkin(4) * t16 + pkin(9) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t60, t46, -t169, -t42, -t72, -t11, -t12, 0, 0;];
tauJ_reg = t5;
