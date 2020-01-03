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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:59:12
% EndTime: 2020-01-03 11:59:14
% DurationCPUTime: 0.67s
% Computational Cost: add. (2661->141), mult. (3756->184), div. (0->0), fcn. (2015->8), ass. (0->100)
t99 = qJD(1) + qJD(2);
t137 = (qJD(5) * t99);
t147 = 2 * t137;
t106 = sin(qJ(4));
t109 = cos(qJ(4));
t97 = t99 ^ 2;
t86 = t109 * t97 * t106;
t80 = qJDD(4) + t86;
t146 = pkin(4) * t80;
t104 = sin(pkin(8));
t105 = cos(pkin(8));
t107 = sin(qJ(2));
t110 = cos(qJ(2));
t108 = sin(qJ(1));
t111 = cos(qJ(1));
t121 = -t111 * g(2) - t108 * g(3);
t77 = qJDD(1) * pkin(1) + t121;
t120 = t108 * g(2) - t111 * g(3);
t78 = -qJD(1) ^ 2 * pkin(1) - t120;
t47 = -t107 * t78 + t110 * t77;
t98 = qJDD(1) + qJDD(2);
t117 = t98 * pkin(2) + t47;
t48 = t107 * t77 + t110 * t78;
t45 = -t97 * pkin(2) + t48;
t34 = t104 * t117 + t105 * t45;
t31 = -t97 * pkin(3) + t98 * pkin(7) + t34;
t143 = t106 * t31;
t102 = -g(1) + qJDD(3);
t94 = t109 * t102;
t24 = -t94 + t143;
t25 = t106 * t102 + t109 * t31;
t9 = t106 * t24 + t109 * t25;
t138 = qJD(4) * t99;
t129 = t106 * t138;
t92 = t109 * t98;
t68 = t92 - t129;
t139 = qJ(5) * t106;
t79 = qJD(4) * pkin(4) - t99 * t139;
t115 = t68 * qJ(5) - qJD(4) * t79 + t109 * t147 + t25;
t100 = t106 ^ 2;
t145 = t100 * t97;
t101 = t109 ^ 2;
t144 = t101 * t97;
t142 = t106 * t80;
t141 = t106 * t98;
t81 = qJDD(4) - t86;
t140 = t109 * t81;
t136 = t100 + t101;
t33 = -t104 * t45 + t105 * t117;
t30 = -t98 * pkin(3) - t97 * pkin(7) - t33;
t4 = t104 * t9 - t105 * t30;
t135 = pkin(2) * t4 - pkin(3) * t30 + pkin(7) * t9;
t112 = qJD(4) ^ 2;
t84 = -t112 - t144;
t56 = t109 * t84 - t142;
t69 = t92 - 0.2e1 * t129;
t37 = t104 * t56 + t105 * t69;
t134 = pkin(2) * t37 + pkin(3) * t69 + pkin(7) * t56;
t83 = -t112 - t145;
t58 = -t106 * t83 - t140;
t128 = t109 * t138;
t66 = 0.2e1 * t128 + t141;
t38 = t104 * t58 - t105 * t66;
t133 = pkin(2) * t38 - pkin(3) * t66 + pkin(7) * t58;
t72 = -t104 * t98 - t105 * t97;
t132 = pkin(2) * t72 - t34;
t74 = t136 * t98;
t75 = t136 * t97;
t44 = t104 * t74 + t105 * t75;
t131 = pkin(2) * t44 + pkin(3) * t75 + pkin(7) * t74;
t127 = t106 * t30 + t133;
t126 = -t109 * t30 + t134;
t119 = t104 * t97 - t105 * t98;
t125 = -pkin(2) * t119 + t33;
t67 = t128 + t141;
t116 = -t94 + (-t128 + t67) * qJ(5) - t146;
t124 = t106 * ((qJ(5) * t98 + t147 + t31) * t106 + t116) + t109 * (qJ(5) * t92 + (t75 - t144) * pkin(4) + t115) + t131;
t19 = t106 * t99 * t79 - t68 * pkin(4) - qJ(5) * t144 + qJDD(5) + t30;
t123 = t106 * (-qJ(5) * t83 + t19) + t109 * (-pkin(4) * t66 - qJ(5) * t81) + t133;
t122 = t131 + t9;
t16 = -0.2e1 * t106 * t137 - t116 - t143;
t17 = -pkin(4) * t144 + t115;
t6 = -t106 * t16 + t109 * t17;
t2 = t104 * t6 - t105 * t19;
t114 = pkin(2) * t2 + pkin(7) * t6 - t16 * t139 - pkin(3) * t19 + t109 * (-pkin(4) * t19 + qJ(5) * t17);
t113 = -t80 * t139 + t109 * (pkin(4) * t69 + qJ(5) * t84 - t19) + t134;
t76 = (t100 - t101) * t97;
t57 = t140 + t106 * (-t112 + t144);
t55 = t109 * (t112 - t145) + t142;
t54 = -t106 * t81 + t109 * t83;
t53 = t106 * t84 + t109 * t80;
t50 = (t68 - t129) * t109;
t49 = (t67 + t128) * t106;
t42 = t106 * t69 + t109 * t66;
t32 = pkin(1) * (t107 * (-t104 * t75 + t105 * t74) + t110 * t44);
t21 = pkin(1) * (t107 * (t104 * t66 + t105 * t58) + t110 * t38);
t20 = pkin(1) * (t107 * (-t104 * t69 + t105 * t56) + t110 * t37);
t12 = t104 * t34 + t105 * t33;
t11 = pkin(2) * t12;
t1 = [0, 0, 0, 0, 0, qJDD(1), t121, t120, 0, 0, 0, 0, 0, 0, 0, t98, pkin(1) * (-t107 * t97 + t110 * t98) + t47, pkin(1) * (-t107 * t98 - t110 * t97) - t48, 0, pkin(1) * (t107 * t48 + t110 * t47), 0, 0, 0, 0, 0, t98, pkin(1) * (t107 * t72 - t110 * t119) + t125, pkin(1) * (t107 * t119 + t110 * t72) + t132, 0, pkin(1) * (t107 * (-t104 * t33 + t105 * t34) + t110 * t12) + t11, t49, t42, t55, t50, t57, 0, t20 + t126, t21 + t127, t32 + t122, pkin(1) * (t107 * (t104 * t30 + t105 * t9) + t110 * t4) + t135, t49, t42, t55, t50, t57, 0, t113 + t20, t21 + t123, t32 + t124, pkin(1) * (t107 * (t104 * t19 + t105 * t6) + t110 * t2) + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t47, -t48, 0, 0, 0, 0, 0, 0, 0, t98, t125, t132, 0, t11, t49, t42, t55, t50, t57, 0, t126, t127, t122, t135, t49, t42, t55, t50, t57, 0, t113, t123, t124, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, 0, t53, t54, 0, t106 * t25 - t109 * t24, 0, 0, 0, 0, 0, 0, t53, t54, 0, t106 * t17 + t109 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t76, t141, t86, t92, qJDD(4), -t24, -t25, 0, 0, -t86, t76, t141, t86, t92, qJDD(4), t16 + t146, (t83 + t144) * pkin(4) - t115, -pkin(4) * t141, pkin(4) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t66, -t75, t19;];
tauJ_reg = t1;
