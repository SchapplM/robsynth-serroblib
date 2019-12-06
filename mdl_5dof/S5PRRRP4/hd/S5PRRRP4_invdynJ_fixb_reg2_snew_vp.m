% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRP4
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:20
% EndTime: 2019-12-05 16:46:23
% DurationCPUTime: 0.74s
% Computational Cost: add. (1779->130), mult. (2457->161), div. (0->0), fcn. (1448->8), ass. (0->96)
t100 = sin(qJ(3));
t103 = cos(qJ(3));
t102 = cos(qJ(4));
t90 = qJD(2) + qJD(3);
t88 = t90 ^ 2;
t99 = sin(qJ(4));
t81 = t99 * t88 * t102;
t74 = qJDD(4) - t81;
t128 = t102 * t74;
t105 = qJD(4) ^ 2;
t93 = t99 ^ 2;
t140 = t93 * t88;
t77 = t105 + t140;
t47 = -t99 * t77 + t128;
t125 = qJD(4) * t90;
t89 = qJDD(2) + qJDD(3);
t135 = t99 * t89;
t61 = 0.2e1 * t102 * t125 + t135;
t32 = t100 * t47 + t103 * t61;
t151 = pkin(2) * t32;
t101 = sin(qJ(2));
t104 = cos(qJ(2));
t150 = t104 * t32 + t101 * (-t100 * t61 + t103 * t47);
t148 = pkin(3) * t61 + pkin(7) * t47;
t146 = qJ(5) * t61;
t97 = sin(pkin(8));
t98 = cos(pkin(8));
t76 = -t98 * g(1) - t97 * g(2);
t95 = -g(3) + qJDD(1);
t52 = -t101 * t76 + t104 * t95;
t110 = qJDD(2) * pkin(2) + t52;
t106 = qJD(2) ^ 2;
t53 = t101 * t95 + t104 * t76;
t50 = -t106 * pkin(2) + t53;
t29 = t100 * t110 + t103 * t50;
t27 = -t88 * pkin(3) + t89 * pkin(7) + t29;
t75 = -t97 * g(1) + t98 * g(2);
t72 = t102 * t75;
t19 = t99 * t27 - t72;
t20 = t102 * t27 + t99 * t75;
t6 = t102 * t20 + t99 * t19;
t94 = t102 ^ 2;
t139 = t94 * t88;
t43 = t128 + t99 * (-t105 + t139);
t131 = qJ(5) * t99;
t114 = -pkin(4) * t102 - t131;
t142 = t114 * t88;
t15 = -qJDD(4) * pkin(4) - t105 * qJ(5) + (t27 + t142) * t99 + qJDD(5) - t72;
t145 = 2 * qJD(5);
t118 = t100 * t50 - t103 * t110;
t26 = -t89 * pkin(3) - t88 * pkin(7) + t118;
t143 = -pkin(3) * t26 + pkin(7) * t6;
t141 = t90 * t99;
t73 = qJDD(4) + t81;
t138 = t99 * t73;
t79 = -t105 - t139;
t46 = t102 * t79 - t138;
t122 = t99 * t125;
t127 = t102 * t89;
t62 = -0.2e1 * t122 + t127;
t134 = pkin(3) * t62 + pkin(7) * t46;
t132 = t93 + t94;
t64 = t132 * t89;
t69 = t132 * t88;
t133 = pkin(3) * t69 + pkin(7) * t64;
t129 = t102 * t61;
t124 = t99 * t26 - t148;
t123 = -t102 * t26 + t134;
t115 = qJDD(4) * qJ(5) + (qJD(4) * t145) + t102 * t142 + t20;
t119 = t99 * (qJ(5) * t69 + t15) + t102 * ((-t105 + t69) * pkin(4) + t115) + t133;
t117 = t133 + t6;
t34 = t99 * t62 + t129;
t67 = -t100 * t89 - t103 * t88;
t113 = t100 * t88 - t103 * t89;
t109 = t26 - (-t122 + t127) * pkin(4) - t146;
t107 = t141 * t145 - t109;
t112 = pkin(4) * t129 + t99 * (-pkin(4) * t122 + t107 + t146) + t148;
t111 = t62 * t131 + t134 + t102 * ((t62 - t122) * pkin(4) + t107);
t12 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t141 + t109;
t14 = -t105 * pkin(4) + t115;
t4 = t102 * t14 + t99 * t15;
t108 = pkin(7) * t4 + (-pkin(3) + t114) * t12;
t70 = (t93 - t94) * t88;
t45 = t138 + t102 * (t105 - t140);
t39 = t61 * t99;
t38 = t62 * t102;
t37 = t100 * t64 + t103 * t69;
t36 = pkin(2) * t37;
t31 = t100 * t46 + t103 * t62;
t30 = pkin(2) * t31;
t21 = t101 * (-t100 * t69 + t103 * t64) + t104 * t37;
t16 = t101 * (-t100 * t62 + t103 * t46) + t104 * t31;
t13 = t100 * t29 - t103 * t118;
t2 = t100 * t6 - t103 * t26;
t1 = t100 * t4 - t103 * t12;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, t104 * qJDD(2) - t101 * t106, -t101 * qJDD(2) - t104 * t106, 0, t101 * t53 + t104 * t52, 0, 0, 0, 0, 0, 0, t101 * t67 - t104 * t113, t101 * t113 + t104 * t67, 0, t101 * (t100 * t118 + t103 * t29) + t104 * t13, 0, 0, 0, 0, 0, 0, t16, -t150, t21, t101 * (t100 * t26 + t103 * t6) + t104 * t2, 0, 0, 0, 0, 0, 0, t16, t21, t150, t101 * (t100 * t12 + t103 * t4) + t104 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t52, -t53, 0, 0, 0, 0, 0, 0, 0, t89, -pkin(2) * t113 - t118, pkin(2) * t67 - t29, 0, pkin(2) * t13, t39, t34, t45, t38, t43, 0, t30 + t123, t124 - t151, t36 + t117, pkin(2) * t2 + t143, t39, t45, -t34, 0, -t43, t38, t111 + t30, t36 + t119, t112 + t151, pkin(2) * t1 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -t118, -t29, 0, 0, t39, t34, t45, t38, t43, 0, t123, t124, t117, t143, t39, t45, -t34, 0, -t43, t38, t111, t119, t112, t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t70, t135, t81, t127, qJDD(4), -t19, -t20, 0, 0, -t81, t135, -t70, qJDD(4), -t127, t81, pkin(4) * t73 + qJ(5) * t79 - t15, (-pkin(4) * t99 + qJ(5) * t102) * t89, qJ(5) * t74 + (-t105 + t77) * pkin(4) + t115, -pkin(4) * t15 + qJ(5) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t135, -t77, t15;];
tauJ_reg = t3;
