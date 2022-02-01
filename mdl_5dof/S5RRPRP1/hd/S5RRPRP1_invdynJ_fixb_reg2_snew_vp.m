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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:20:00
% EndTime: 2022-01-20 10:20:03
% DurationCPUTime: 0.72s
% Computational Cost: add. (2661->141), mult. (3756->184), div. (0->0), fcn. (2015->8), ass. (0->101)
t100 = qJD(1) + qJD(2);
t137 = (qJD(5) * t100);
t149 = 2 * t137;
t107 = sin(qJ(4));
t110 = cos(qJ(4));
t98 = t100 ^ 2;
t86 = t110 * t98 * t107;
t80 = qJDD(4) + t86;
t148 = pkin(4) * t80;
t105 = sin(pkin(8));
t106 = cos(pkin(8));
t108 = sin(qJ(2));
t111 = cos(qJ(2));
t109 = sin(qJ(1));
t112 = cos(qJ(1));
t128 = t109 * g(1) - t112 * g(2);
t77 = qJDD(1) * pkin(1) + t128;
t121 = t112 * g(1) + t109 * g(2);
t78 = -qJD(1) ^ 2 * pkin(1) - t121;
t47 = -t108 * t78 + t111 * t77;
t99 = qJDD(1) + qJDD(2);
t118 = t99 * pkin(2) + t47;
t48 = t108 * t77 + t111 * t78;
t45 = -t98 * pkin(2) + t48;
t34 = t105 * t118 + t106 * t45;
t31 = -t98 * pkin(3) + t99 * pkin(7) + t34;
t145 = t107 * t31;
t103 = -g(3) + qJDD(3);
t94 = t110 * t103;
t24 = -t94 + t145;
t25 = t107 * t103 + t110 * t31;
t9 = t107 * t24 + t110 * t25;
t138 = qJD(4) * t100;
t130 = t107 * t138;
t92 = t110 * t99;
t68 = t92 - t130;
t140 = t100 * t107;
t79 = qJD(4) * pkin(4) - qJ(5) * t140;
t116 = t68 * qJ(5) - qJD(4) * t79 + t110 * t149 + t25;
t101 = t107 ^ 2;
t147 = t101 * t98;
t102 = t110 ^ 2;
t146 = t102 * t98;
t144 = t107 * t80;
t143 = t107 * t99;
t81 = qJDD(4) - t86;
t142 = t110 * t81;
t141 = qJ(5) * t107;
t139 = t101 + t102;
t33 = -t105 * t45 + t106 * t118;
t30 = -t99 * pkin(3) - t98 * pkin(7) - t33;
t4 = t105 * t9 - t106 * t30;
t136 = pkin(2) * t4 - pkin(3) * t30 + pkin(7) * t9;
t113 = qJD(4) ^ 2;
t84 = -t113 - t146;
t57 = t110 * t84 - t144;
t69 = t92 - 0.2e1 * t130;
t37 = t105 * t57 + t106 * t69;
t135 = pkin(2) * t37 + pkin(3) * t69 + pkin(7) * t57;
t83 = -t113 - t147;
t58 = -t107 * t83 - t142;
t129 = t110 * t138;
t66 = 0.2e1 * t129 + t143;
t38 = t105 * t58 - t106 * t66;
t134 = pkin(2) * t38 - pkin(3) * t66 + pkin(7) * t58;
t72 = -t105 * t99 - t106 * t98;
t133 = pkin(2) * t72 - t34;
t74 = t139 * t99;
t75 = t139 * t98;
t44 = t105 * t74 + t106 * t75;
t132 = pkin(2) * t44 + pkin(3) * t75 + pkin(7) * t74;
t127 = t107 * t30 + t134;
t126 = -t110 * t30 + t135;
t120 = t105 * t98 - t106 * t99;
t125 = -pkin(2) * t120 + t33;
t67 = t129 + t143;
t117 = -t94 + (-t129 + t67) * qJ(5) - t148;
t124 = t107 * ((qJ(5) * t99 + t149 + t31) * t107 + t117) + t110 * (qJ(5) * t92 + (t75 - t146) * pkin(4) + t116) + t132;
t19 = -t68 * pkin(4) - qJ(5) * t146 + t79 * t140 + qJDD(5) + t30;
t123 = t107 * (-qJ(5) * t83 + t19) + t110 * (-pkin(4) * t66 - qJ(5) * t81) + t134;
t122 = t132 + t9;
t16 = -0.2e1 * t107 * t137 - t117 - t145;
t17 = -pkin(4) * t146 + t116;
t6 = -t107 * t16 + t110 * t17;
t2 = t105 * t6 - t106 * t19;
t115 = pkin(2) * t2 + pkin(7) * t6 - t16 * t141 - pkin(3) * t19 + t110 * (-pkin(4) * t19 + qJ(5) * t17);
t114 = -t80 * t141 + t110 * (pkin(4) * t69 + qJ(5) * t84 - t19) + t135;
t76 = (t101 - t102) * t98;
t56 = -t107 * t81 + t110 * t83;
t55 = t144 + t110 * (t113 - t147);
t54 = t107 * t84 + t110 * t80;
t53 = t107 * (-t113 + t146) + t142;
t50 = (t67 + t129) * t107;
t49 = (t68 - t130) * t110;
t42 = t107 * t69 + t110 * t66;
t32 = pkin(1) * (t108 * (-t105 * t75 + t106 * t74) + t111 * t44);
t21 = pkin(1) * (t108 * (t105 * t66 + t106 * t58) + t111 * t38);
t20 = pkin(1) * (t108 * (-t105 * t69 + t106 * t57) + t111 * t37);
t12 = t105 * t34 + t106 * t33;
t11 = pkin(2) * t12;
t1 = [0, 0, 0, 0, 0, qJDD(1), t128, t121, 0, 0, 0, 0, 0, 0, 0, t99, pkin(1) * (-t108 * t98 + t111 * t99) + t47, pkin(1) * (-t108 * t99 - t111 * t98) - t48, 0, pkin(1) * (t108 * t48 + t111 * t47), 0, 0, 0, 0, 0, t99, pkin(1) * (t108 * t72 - t111 * t120) + t125, pkin(1) * (t108 * t120 + t111 * t72) + t133, 0, pkin(1) * (t108 * (-t105 * t33 + t106 * t34) + t111 * t12) + t11, t50, t42, t55, t49, t53, 0, t20 + t126, t21 + t127, t32 + t122, pkin(1) * (t108 * (t105 * t30 + t106 * t9) + t111 * t4) + t136, t50, t42, t55, t49, t53, 0, t114 + t20, t21 + t123, t32 + t124, pkin(1) * (t108 * (t105 * t19 + t106 * t6) + t111 * t2) + t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t47, -t48, 0, 0, 0, 0, 0, 0, 0, t99, t125, t133, 0, t11, t50, t42, t55, t49, t53, 0, t126, t127, t122, t136, t50, t42, t55, t49, t53, 0, t114, t123, t124, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, 0, 0, 0, 0, 0, t54, t56, 0, t107 * t25 - t110 * t24, 0, 0, 0, 0, 0, 0, t54, t56, 0, t107 * t17 + t110 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t76, t143, t86, t92, qJDD(4), -t24, -t25, 0, 0, -t86, t76, t143, t86, t92, qJDD(4), t16 + t148, (t83 + t146) * pkin(4) - t116, -pkin(4) * t143, pkin(4) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t66, -t75, t19;];
tauJ_reg = t1;
