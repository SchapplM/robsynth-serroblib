% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:14
% EndTime: 2019-12-05 16:44:18
% DurationCPUTime: 1.36s
% Computational Cost: add. (3247->203), mult. (6822->261), div. (0->0), fcn. (4597->8), ass. (0->134)
t117 = qJDD(3) + qJDD(4);
t126 = sin(qJ(4));
t129 = cos(qJ(4));
t130 = cos(qJ(3));
t127 = sin(qJ(3));
t151 = qJD(2) * t127;
t91 = -t129 * t130 * qJD(2) + t126 * t151;
t93 = (t130 * t126 + t127 * t129) * qJD(2);
t75 = t93 * t91;
t177 = t75 - t117;
t178 = t177 * pkin(4);
t176 = t126 * t177;
t175 = t129 * t177;
t132 = qJD(2) ^ 2;
t107 = t127 * t132 * t130;
t102 = qJDD(3) + t107;
t149 = qJD(2) * qJD(3);
t144 = t130 * t149;
t121 = -g(3) + qJDD(1);
t123 = sin(pkin(8));
t124 = cos(pkin(8));
t101 = -t124 * g(1) - t123 * g(2);
t128 = sin(qJ(2));
t139 = t123 * g(1) - t124 * g(2);
t168 = cos(qJ(2));
t133 = -t168 * t101 - t128 * t139;
t156 = qJDD(2) * pkin(6);
t72 = -t132 * pkin(2) - t133 + t156;
t64 = -t130 * t121 + t127 * t72;
t112 = t127 * qJDD(2);
t98 = t112 + t144;
t41 = (-t98 + t144) * pkin(7) + t102 * pkin(3) - t64;
t104 = qJD(3) * pkin(3) - pkin(7) * t151;
t120 = t130 ^ 2;
t115 = t120 * t132;
t65 = t127 * t121 + t130 * t72;
t113 = t130 * qJDD(2);
t145 = t127 * t149;
t99 = t113 - t145;
t42 = -pkin(3) * t115 + t99 * pkin(7) - qJD(3) * t104 + t65;
t20 = t126 * t42 - t129 * t41;
t118 = qJD(3) + qJD(4);
t83 = t118 * t91;
t173 = qJ(5) * t83 + 0.2e1 * qJD(5) * t93 + t178 + t20;
t23 = t126 * t41 + t129 * t42;
t143 = t126 * t98 - t129 * t99;
t60 = -t93 * qJD(4) - t143;
t79 = t118 * pkin(4) - t93 * qJ(5);
t89 = t91 ^ 2;
t14 = -t89 * pkin(4) + t60 * qJ(5) - 0.2e1 * qJD(5) * t91 - t118 * t79 + t23;
t90 = t93 ^ 2;
t116 = t118 ^ 2;
t135 = (-qJD(4) + t118) * t93 - t143;
t138 = t126 * t99 + t129 * t98;
t61 = -t91 * qJD(4) + t138;
t53 = t61 + t83;
t27 = t126 * t135 - t129 * t53;
t172 = pkin(7) * t27;
t67 = -t116 - t89;
t36 = t126 * t67 - t175;
t171 = pkin(7) * t36;
t69 = t75 + t117;
t163 = t126 * t69;
t78 = -t90 - t116;
t54 = t129 * t78 - t163;
t170 = pkin(7) * t54;
t28 = t126 * t53 + t129 * t135;
t62 = -t89 - t90;
t169 = pkin(6) * (-t127 * t27 + t130 * t28) - pkin(2) * t62;
t6 = t126 * t23 - t129 * t20;
t167 = t127 * t6;
t37 = t129 * t67 + t176;
t150 = qJD(4) + t118;
t48 = t150 * t93 + t143;
t166 = pkin(6) * (-t127 * t36 + t130 * t37) - pkin(2) * t48;
t51 = -t150 * t91 + t138;
t161 = t129 * t69;
t55 = -t126 * t78 - t161;
t165 = pkin(6) * (-t127 * t54 + t130 * t55) - pkin(2) * t51;
t141 = -t128 * t101 + t168 * t139;
t71 = -qJDD(2) * pkin(2) - t132 * pkin(6) - t141;
t56 = -t99 * pkin(3) - pkin(7) * t115 + t104 * t151 + t71;
t164 = t126 * t56;
t162 = t129 * t56;
t160 = qJ(5) * t126;
t159 = qJ(5) * t129;
t155 = t118 * t126;
t154 = t118 * t129;
t153 = t127 * t102;
t103 = qJDD(3) - t107;
t152 = t130 * t103;
t148 = -pkin(3) * t62 + pkin(7) * t28;
t147 = -pkin(3) * t48 + pkin(7) * t37;
t146 = -pkin(3) * t51 + pkin(7) * t55;
t7 = t126 * t20 + t129 * t23;
t142 = t127 * t64 + t130 * t65;
t137 = pkin(4) * t78 - t14;
t12 = -t61 * qJ(5) - t173;
t134 = t12 - t178;
t16 = -t60 * pkin(4) - t89 * qJ(5) + t93 * t79 + qJDD(5) + t56;
t131 = qJD(3) ^ 2;
t119 = t127 ^ 2;
t114 = t119 * t132;
t106 = -t115 - t131;
t105 = -t114 - t131;
t100 = t113 - 0.2e1 * t145;
t97 = t112 + 0.2e1 * t144;
t81 = -t90 + t116;
t80 = t89 - t116;
t73 = t90 - t89;
t52 = t61 - t83;
t47 = pkin(3) * t54;
t43 = pkin(4) * t53;
t35 = pkin(3) * t36;
t33 = (t127 * (t126 * t93 - t129 * t91) + t130 * (-t126 * t91 - t129 * t93)) * t118;
t32 = -pkin(4) * t51 - qJ(5) * t69;
t31 = t127 * (t129 * t80 - t163) + t130 * (t126 * t80 + t161);
t30 = t127 * (-t126 * t81 - t175) + t130 * (t129 * t81 - t176);
t29 = t127 * t55 + t130 * t54;
t25 = pkin(3) * t27;
t22 = t127 * (t129 * t61 - t93 * t155) + t130 * (t126 * t61 + t93 * t154);
t21 = t127 * (-t126 * t60 + t91 * t154) + t130 * (t129 * t60 + t91 * t155);
t18 = t127 * t37 + t130 * t36;
t15 = -qJ(5) * t78 + t16;
t13 = -pkin(4) * t48 + qJ(5) * t67 - t16;
t11 = pkin(4) * t12;
t10 = t127 * t28 + t130 * t27;
t9 = t127 * (-t126 * t52 - t129 * t48) + t130 * (-t126 * t48 + t129 * t52);
t5 = (t53 + t61) * qJ(5) + t173;
t4 = -pkin(4) * t62 + qJ(5) * t135 + t14;
t3 = -pkin(4) * t16 + qJ(5) * t14;
t2 = -t126 * t12 + t129 * t14;
t1 = t129 * t12 + t126 * t14;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, 0, 0, 0, 0, 0, t130 * t102 + t127 * t106, -t127 * t103 + t130 * t105, 0, t127 * t65 - t130 * t64, 0, 0, 0, 0, 0, 0, t18, t29, t10, t127 * t7 + t130 * t6, 0, 0, 0, 0, 0, 0, t18, t29, t10, t130 * t1 + t127 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t141, t133, 0, 0, (t98 + t144) * t127, t127 * t100 + t130 * t97, t153 + t130 * (-t114 + t131), (t99 - t145) * t130, t127 * (t115 - t131) + t152, 0, -t130 * t71 + pkin(2) * t100 + pkin(6) * (t106 * t130 - t153), t127 * t71 - pkin(2) * t97 + pkin(6) * (-t105 * t127 - t152), pkin(2) * (t114 + t115) + (t119 + t120) * t156 + t142, -pkin(2) * t71 + pkin(6) * t142, t22, t9, t30, t21, t31, t33, t127 * (t164 - t171) + t130 * (t147 - t162) + t166, t127 * (t162 - t170) + t130 * (t146 + t164) + t165, t127 * (-t6 - t172) + t130 * (t148 + t7) + t169, -pkin(7) * t167 + t130 * (-pkin(3) * t56 + pkin(7) * t7) - pkin(2) * t56 + pkin(6) * (t130 * t7 - t167), t22, t9, t30, t21, t31, t33, t127 * (-t126 * t13 + t159 * t177 - t171) + t130 * (t129 * t13 + t160 * t177 + t147) + t166, t127 * (-t126 * t32 + t129 * t15 - t170) + t130 * (t126 * t15 + t129 * t32 + t146) + t165, t127 * (-t126 * t4 + t129 * t5 - t172) + t130 * (t126 * t5 + t129 * t4 + t148) + t169, t127 * (-pkin(7) * t1 - t12 * t159 - t126 * t3) + t130 * (-pkin(3) * t16 + pkin(7) * t2 - t12 * t160 + t129 * t3) - pkin(2) * t16 + pkin(6) * (-t1 * t127 + t130 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, t114 - t115, t112, t107, t113, qJDD(3), -t64, -t65, 0, 0, t75, t73, t53, -t75, t135, t117, -t20 + t35, t47 - t23, t25, pkin(3) * t6, t75, t73, t53, -t75, t135, t117, t134 + t35, t47 + t137, -t43 + t25, pkin(3) * t1 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t73, t53, -t75, t135, t117, -t20, -t23, 0, 0, t75, t73, t53, -t75, t135, t117, t134, t137, -t43, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t51, t62, t16;];
tauJ_reg = t8;
