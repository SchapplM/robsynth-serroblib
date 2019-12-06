% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:23
% EndTime: 2019-12-05 18:22:27
% DurationCPUTime: 0.71s
% Computational Cost: add. (2661->141), mult. (3756->184), div. (0->0), fcn. (2015->8), ass. (0->101)
t101 = qJD(1) + qJD(2);
t136 = (qJD(5) * t101);
t149 = 2 * t136;
t108 = sin(qJ(4));
t111 = cos(qJ(4));
t99 = t101 ^ 2;
t86 = t111 * t99 * t108;
t80 = qJDD(4) + t86;
t148 = pkin(4) * t80;
t100 = qJDD(1) + qJDD(2);
t106 = sin(pkin(8));
t107 = cos(pkin(8));
t109 = sin(qJ(2));
t112 = cos(qJ(2));
t110 = sin(qJ(1));
t113 = cos(qJ(1));
t147 = t113 * g(2) + t110 * g(3);
t77 = qJDD(1) * pkin(1) + t147;
t121 = -t110 * g(2) + t113 * g(3);
t78 = -qJD(1) ^ 2 * pkin(1) - t121;
t47 = -t109 * t78 + t112 * t77;
t119 = t100 * pkin(2) + t47;
t48 = t109 * t77 + t112 * t78;
t45 = -t99 * pkin(2) + t48;
t34 = t106 * t119 + t107 * t45;
t31 = -t99 * pkin(3) + t100 * pkin(7) + t34;
t144 = t108 * t31;
t104 = -g(1) + qJDD(3);
t94 = t111 * t104;
t24 = -t94 + t144;
t25 = t108 * t104 + t111 * t31;
t9 = t108 * t24 + t111 * t25;
t137 = qJD(4) * t101;
t129 = t108 * t137;
t92 = t111 * t100;
t68 = t92 - t129;
t140 = t101 * t108;
t79 = qJD(4) * pkin(4) - qJ(5) * t140;
t117 = t68 * qJ(5) - qJD(4) * t79 + t111 * t149 + t25;
t102 = t108 ^ 2;
t146 = t102 * t99;
t103 = t111 ^ 2;
t145 = t103 * t99;
t143 = t108 * t80;
t81 = qJDD(4) - t86;
t142 = t111 * t81;
t141 = qJ(5) * t108;
t139 = t108 * t100;
t138 = t102 + t103;
t33 = -t106 * t45 + t107 * t119;
t30 = -t100 * pkin(3) - t99 * pkin(7) - t33;
t4 = t106 * t9 - t107 * t30;
t135 = pkin(2) * t4 - pkin(3) * t30 + pkin(7) * t9;
t114 = qJD(4) ^ 2;
t84 = -t114 - t145;
t57 = t111 * t84 - t143;
t69 = t92 - 0.2e1 * t129;
t37 = t106 * t57 + t107 * t69;
t134 = pkin(2) * t37 + pkin(3) * t69 + pkin(7) * t57;
t83 = -t114 - t146;
t58 = -t108 * t83 - t142;
t128 = t111 * t137;
t66 = 0.2e1 * t128 + t139;
t38 = t106 * t58 - t107 * t66;
t133 = pkin(2) * t38 - pkin(3) * t66 + pkin(7) * t58;
t72 = -t106 * t100 - t107 * t99;
t132 = pkin(2) * t72 - t34;
t74 = t138 * t100;
t75 = t138 * t99;
t44 = t106 * t74 + t107 * t75;
t131 = pkin(2) * t44 + pkin(3) * t75 + pkin(7) * t74;
t127 = t108 * t30 + t133;
t126 = -t111 * t30 + t134;
t73 = t107 * t100 - t106 * t99;
t125 = pkin(2) * t73 + t33;
t67 = t128 + t139;
t118 = -t94 + (-t128 + t67) * qJ(5) - t148;
t124 = t108 * ((qJ(5) * t100 + t149 + t31) * t108 + t118) + t111 * (qJ(5) * t92 + (t75 - t145) * pkin(4) + t117) + t131;
t19 = -t68 * pkin(4) - qJ(5) * t145 + t140 * t79 + qJDD(5) + t30;
t123 = t108 * (-qJ(5) * t83 + t19) + t111 * (-pkin(4) * t66 - qJ(5) * t81) + t133;
t122 = t131 + t9;
t16 = -0.2e1 * t108 * t136 - t118 - t144;
t17 = -pkin(4) * t145 + t117;
t6 = -t108 * t16 + t111 * t17;
t2 = t106 * t6 - t107 * t19;
t116 = pkin(2) * t2 + pkin(7) * t6 - t141 * t16 - pkin(3) * t19 + t111 * (-pkin(4) * t19 + qJ(5) * t17);
t115 = -t141 * t80 + t111 * (pkin(4) * t69 + qJ(5) * t84 - t19) + t134;
t76 = (t102 - t103) * t99;
t56 = -t108 * t81 + t111 * t83;
t55 = t143 + t111 * (t114 - t146);
t54 = t108 * t84 + t111 * t80;
t53 = t108 * (-t114 + t145) + t142;
t50 = (t67 + t128) * t108;
t49 = (t68 - t129) * t111;
t42 = t108 * t69 + t111 * t66;
t32 = pkin(1) * (t109 * (-t106 * t75 + t107 * t74) + t112 * t44);
t21 = pkin(1) * (t109 * (t106 * t66 + t107 * t58) + t112 * t38);
t20 = pkin(1) * (t109 * (-t106 * t69 + t107 * t57) + t112 * t37);
t12 = t106 * t34 + t107 * t33;
t11 = pkin(2) * t12;
t1 = [0, 0, 0, 0, 0, qJDD(1), t147, t121, 0, 0, 0, 0, 0, 0, 0, t100, pkin(1) * (t112 * t100 - t109 * t99) + t47, pkin(1) * (-t109 * t100 - t112 * t99) - t48, 0, pkin(1) * (t109 * t48 + t112 * t47), 0, 0, 0, 0, 0, t100, pkin(1) * (t109 * t72 + t112 * t73) + t125, pkin(1) * (-t109 * t73 + t112 * t72) + t132, 0, pkin(1) * (t109 * (-t106 * t33 + t107 * t34) + t112 * t12) + t11, t50, t42, t55, t49, t53, 0, t20 + t126, t21 + t127, t32 + t122, pkin(1) * (t109 * (t106 * t30 + t107 * t9) + t112 * t4) + t135, t50, t42, t55, t49, t53, 0, t115 + t20, t21 + t123, t32 + t124, pkin(1) * (t109 * (t106 * t19 + t107 * t6) + t112 * t2) + t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t47, -t48, 0, 0, 0, 0, 0, 0, 0, t100, t125, t132, 0, t11, t50, t42, t55, t49, t53, 0, t126, t127, t122, t135, t50, t42, t55, t49, t53, 0, t115, t123, t124, t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0, 0, 0, 0, 0, t54, t56, 0, t108 * t25 - t111 * t24, 0, 0, 0, 0, 0, 0, t54, t56, 0, t108 * t17 + t111 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t76, t139, t86, t92, qJDD(4), -t24, -t25, 0, 0, -t86, t76, t139, t86, t92, qJDD(4), t16 + t148, (t83 + t145) * pkin(4) - t117, -pkin(4) * t139, pkin(4) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t66, -t75, t19;];
tauJ_reg = t1;
