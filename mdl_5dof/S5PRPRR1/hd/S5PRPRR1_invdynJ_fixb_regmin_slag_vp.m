% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRR1
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:07
% EndTime: 2019-12-05 15:43:11
% DurationCPUTime: 0.85s
% Computational Cost: add. (1068->172), mult. (2362->226), div. (0->0), fcn. (1914->12), ass. (0->110)
t100 = cos(qJ(5));
t101 = cos(qJ(4));
t97 = cos(pkin(9));
t139 = t101 * t97;
t96 = sin(pkin(9));
t99 = sin(qJ(4));
t144 = t99 * t96;
t111 = t139 - t144;
t47 = t111 * qJD(2);
t56 = t101 * t96 + t99 * t97;
t48 = t56 * qJD(2);
t98 = sin(qJ(5));
t19 = -t100 * t47 + t98 * t48;
t95 = qJD(4) + qJD(5);
t146 = t19 * t95;
t133 = qJD(5) * t100;
t137 = qJD(5) * t98;
t126 = qJD(2) * t144;
t128 = qJDD(2) * t101;
t131 = t97 * qJDD(2);
t134 = qJD(4) * t101;
t138 = qJD(2) * t97;
t127 = t96 * t128 + t99 * t131 + t134 * t138;
t26 = -qJD(4) * t126 + t127;
t132 = t96 * qJDD(2);
t120 = -t97 * t128 + t99 * t132;
t50 = t56 * qJD(4);
t27 = qJD(2) * t50 + t120;
t5 = t100 * t26 + t47 * t133 - t48 * t137 - t98 * t27;
t159 = t5 + t146;
t114 = t100 * t48 + t98 * t47;
t158 = t114 * t19;
t147 = t114 * t95;
t6 = t114 * qJD(5) + t100 * t27 + t98 * t26;
t157 = -t6 + t147;
t129 = qJD(2) * qJD(3);
t152 = qJ(3) * qJDD(2) + t129;
t94 = pkin(8) + qJ(2);
t86 = sin(t94);
t88 = cos(t94);
t122 = g(1) * t88 + g(2) * t86;
t156 = t114 ^ 2 - t19 ^ 2;
t143 = pkin(6) + qJ(3);
t65 = t143 * t96;
t84 = t97 * qJD(1);
t44 = -qJD(2) * t65 + t84;
t135 = qJ(3) * qJD(2);
t62 = t96 * qJD(1) + t97 * t135;
t45 = pkin(6) * t138 + t62;
t112 = -t101 * t45 - t99 * t44;
t14 = t47 * pkin(7) - t112;
t82 = t97 * qJDD(1);
t34 = t82 + (-t143 * qJDD(2) - t129) * t96;
t43 = t96 * qJDD(1) + t152 * t97;
t35 = pkin(6) * t131 + t43;
t124 = t101 * t34 - t99 * t35;
t2 = qJDD(4) * pkin(4) - t26 * pkin(7) + t112 * qJD(4) + t124;
t79 = -t97 * pkin(3) - pkin(2);
t64 = t79 * qJD(2) + qJD(3);
t32 = -t47 * pkin(4) + t64;
t93 = pkin(9) + qJ(4);
t89 = qJ(5) + t93;
t77 = sin(t89);
t78 = cos(t89);
t155 = t32 * t19 + g(3) * t77 + t14 * t137 + (-t14 * t95 - t2) * t98 + t122 * t78;
t66 = t143 * t97;
t142 = t101 * t66 - t99 * t65;
t153 = t101 * t44 - t99 * t45;
t13 = -t48 * pkin(7) + t153;
t10 = qJD(4) * pkin(4) + t13;
t140 = t100 * t14;
t117 = -t98 * t10 - t140;
t113 = t101 * t35 + t99 * t34;
t3 = -t27 * pkin(7) + t153 * qJD(4) + t113;
t151 = -g(3) * t78 + t117 * qJD(5) + t100 * t2 - t32 * t114 + t122 * t77 - t98 * t3;
t150 = pkin(4) * t50;
t141 = t96 ^ 2 + t97 ^ 2;
t136 = qJDD(2) * pkin(2);
t123 = -t101 * t65 - t99 * t66;
t121 = g(1) * t86 - g(2) * t88;
t31 = t100 * t56 + t111 * t98;
t30 = -t100 * t111 + t98 * t56;
t49 = t111 * qJD(4);
t7 = -t30 * qJD(5) + t100 * t49 - t98 * t50;
t90 = qJDD(4) + qJDD(5);
t119 = t31 * t90 + t7 * t95;
t118 = (-t96 * t135 + t84) * t96 - t62 * t97;
t16 = -t56 * pkin(7) + t123;
t17 = pkin(7) * t111 + t142;
t116 = t100 * t17 + t98 * t16;
t115 = t100 * t16 - t98 * t17;
t63 = t79 * qJDD(2) + qJDD(3);
t110 = -t121 - t136;
t80 = qJDD(3) - t136;
t109 = -t110 - t80;
t107 = -t65 * t134 + qJD(3) * t139 + (-qJD(3) * t96 - qJD(4) * t66) * t99;
t42 = -t152 * t96 + t82;
t106 = -t42 * t96 + t43 * t97 - t122;
t103 = -t56 * qJD(3) - t142 * qJD(4);
t87 = cos(t93);
t85 = sin(t93);
t37 = -pkin(4) * t111 + t79;
t29 = -t50 * qJD(4) + qJDD(4) * t111;
t28 = t49 * qJD(4) + t56 * qJDD(4);
t15 = t27 * pkin(4) + t63;
t12 = -t49 * pkin(7) + t103;
t11 = -t50 * pkin(7) + t107;
t8 = t31 * qJD(5) + t100 * t50 + t98 * t49;
t4 = -t30 * t90 - t8 * t95;
t1 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, t42 * t97 + t43 * t96 - g(3), 0, 0, 0, 0, 0, t29, -t28, 0, 0, 0, 0, 0, t4, -t119; 0, qJDD(2), t121, t122, t109 * t97, -t109 * t96, t152 * t141 + t106, -t118 * qJD(3) + (t121 - t80) * pkin(2) + t106 * qJ(3), t26 * t56 + t48 * t49, t111 * t26 - t56 * t27 + t49 * t47 - t48 * t50, t28, t29, 0, t103 * qJD(4) + t123 * qJDD(4) - t111 * t63 + t121 * t87 + t79 * t27 + t64 * t50, -t107 * qJD(4) - t142 * qJDD(4) - t121 * t85 + t79 * t26 + t64 * t49 + t63 * t56, t114 * t7 + t5 * t31, -t114 * t8 - t7 * t19 - t5 * t30 - t31 * t6, t119, t4, 0, t19 * t150 + t37 * t6 + t15 * t30 + t32 * t8 + (-t116 * qJD(5) + t100 * t12 - t98 * t11) * t95 + t115 * t90 + t121 * t78, t114 * t150 + t37 * t5 + t15 * t31 + t32 * t7 - (t115 * qJD(5) + t100 * t11 + t98 * t12) * t95 - t116 * t90 - t121 * t77; 0, 0, 0, 0, -t131, t132, -t141 * qJD(2) ^ 2, t118 * qJD(2) + qJDD(3) + t110, 0, 0, 0, 0, 0, 0.2e1 * t48 * qJD(4) + t120, (t47 - t126) * qJD(4) + t127, 0, 0, 0, 0, 0, t6 + t147, t5 - t146; 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t47, -t47 ^ 2 + t48 ^ 2, (-t47 - t126) * qJD(4) + t127, -t120, qJDD(4), -g(3) * t87 + t122 * t85 - t64 * t48 + t124, g(3) * t85 + t122 * t87 - t64 * t47 - t113, t158, t156, t159, t157, t90, -(-t98 * t13 - t140) * t95 + (t100 * t90 - t95 * t137 - t48 * t19) * pkin(4) + t151, (-qJD(5) * t10 + t13 * t95 - t3) * t100 + (-t114 * t48 - t95 * t133 - t98 * t90) * pkin(4) + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t156, t159, t157, t90, -t117 * t95 + t151, (-t3 + (-qJD(5) + t95) * t10) * t100 + t155;];
tau_reg = t1;
