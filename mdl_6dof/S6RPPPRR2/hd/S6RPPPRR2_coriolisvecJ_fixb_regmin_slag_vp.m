% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:02
% EndTime: 2019-03-09 01:32:05
% DurationCPUTime: 0.89s
% Computational Cost: add. (1198->163), mult. (2659->233), div. (0->0), fcn. (1930->8), ass. (0->95)
t132 = 2 * qJD(3);
t73 = sin(pkin(10));
t75 = cos(pkin(10));
t78 = sin(qJ(5));
t80 = cos(qJ(5));
t50 = t80 * t73 + t78 * t75;
t86 = qJD(1) * t50;
t129 = qJD(6) + t86;
t79 = cos(qJ(6));
t103 = t79 * qJD(5);
t108 = qJD(1) * t75;
t98 = t80 * t108;
t109 = qJD(1) * t73;
t99 = t78 * t109;
t45 = t98 - t99;
t77 = sin(qJ(6));
t26 = t77 * t45 - t103;
t131 = t129 * t26;
t60 = sin(pkin(9)) * pkin(1) + qJ(3);
t130 = qJD(1) * t60;
t54 = qJD(5) * t99;
t34 = qJD(5) * t98 - t54;
t112 = t77 * t34;
t94 = t79 * t129;
t128 = -t129 * t94 - t112;
t126 = t86 * qJD(5);
t59 = -cos(pkin(9)) * pkin(1) - pkin(2) - qJ(4);
t125 = t59 * qJD(1);
t124 = -qJD(6) + t129;
t48 = qJD(3) + t125;
t29 = -t73 * qJD(2) + t75 * t48;
t24 = -pkin(7) * t108 + t29;
t30 = t75 * qJD(2) + t73 * t48;
t25 = -pkin(7) * t109 + t30;
t8 = t78 * t24 + t80 * t25;
t49 = t78 * t73 - t80 * t75;
t85 = t49 * qJD(4);
t4 = -qJD(1) * t85 + t8 * qJD(5);
t123 = (t45 * pkin(5) + t129 * pkin(8)) * t129 + t4;
t32 = t79 * t34;
t104 = qJD(6) * t77;
t106 = qJD(5) * t80;
t107 = qJD(5) * t78;
t46 = -t73 * t106 - t75 * t107;
t88 = t49 * t104 + t79 * t46;
t122 = -t129 * t88 + t49 * t32;
t120 = -pkin(7) + t59;
t41 = t120 * t73;
t42 = t120 * t75;
t19 = t78 * t41 - t80 * t42;
t87 = t50 * qJD(4);
t10 = -t19 * qJD(5) - t87;
t51 = qJD(4) + t130;
t40 = pkin(4) * t109 + t51;
t12 = pkin(5) * t86 - t45 * pkin(8) + t40;
t52 = t73 * pkin(4) + t60;
t18 = t50 * pkin(5) + t49 * pkin(8) + t52;
t20 = t80 * t41 + t78 * t42;
t7 = t80 * t24 - t78 * t25;
t3 = -qJD(1) * t87 + t7 * qJD(5);
t5 = -qJD(5) * pkin(5) - t7;
t121 = -(qJD(6) * t18 + t10) * t129 - t20 * t34 - (qJD(6) * t12 + t3) * t50 - t4 * t49 + t5 * t46;
t71 = qJD(3) * qJD(1);
t97 = 0.2e1 * t71;
t13 = qJD(6) * t103 - t45 * t104 - t126 * t79;
t28 = t77 * qJD(5) + t79 * t45;
t47 = t75 * t106 - t73 * t107;
t119 = t13 * t50 + t28 * t47;
t118 = t13 * t77;
t117 = t18 * t34;
t116 = t28 * t45;
t115 = t45 * t26;
t114 = t49 * t13;
t113 = t77 * t126;
t111 = t73 ^ 2 + t75 ^ 2;
t105 = qJD(6) * t49;
t36 = t47 * qJD(5);
t101 = qJD(1) * qJD(4);
t92 = qJD(6) * t50 + qJD(1);
t6 = qJD(5) * pkin(8) + t8;
t1 = t79 * t12 - t77 * t6;
t2 = t77 * t12 + t79 * t6;
t14 = t28 * qJD(6) - t113;
t91 = -t50 * t14 - t47 * t26;
t90 = t29 * t75 + t30 * t73;
t89 = t32 + (-t77 * t86 - t104) * t129;
t83 = -pkin(8) * t34 + (t5 + t7) * t129;
t82 = (t79 * t105 - t77 * t46) * t129 + t49 * t112;
t81 = qJD(1) ^ 2;
t35 = t46 * qJD(5);
t21 = t47 * pkin(5) - t46 * pkin(8) + qJD(3);
t16 = t34 * pkin(5) + pkin(8) * t126 + t71;
t15 = t79 * t16;
t11 = t20 * qJD(5) - t85;
t9 = [0, 0, 0, 0, 0, t97, t130 * t132, t73 * t97, t75 * t97, 0.2e1 * t111 * t101 (t51 + t130) * qJD(3) + (-t111 * t125 - t90) * qJD(4), t126 * t49 + t45 * t46, t126 * t50 + t49 * t34 - t45 * t47 - t46 * t86, t35, -t36, 0, -t11 * qJD(5) + t86 * t132 + t52 * t34 + t40 * t47, -t10 * qJD(5) - t52 * t126 + t40 * t46 + (-qJD(1) * t49 + t45) * qJD(3), -t79 * t114 + t88 * t28 (-t26 * t79 - t28 * t77) * t46 + (t118 + t14 * t79 + (-t26 * t77 + t28 * t79) * qJD(6)) * t49, t119 - t122, t82 + t91, t129 * t47 + t34 * t50, t1 * t47 + t11 * t26 + t19 * t14 + t15 * t50 + (t117 + t21 * t129 + (-t129 * t20 - t5 * t49 - t6 * t50) * qJD(6)) * t79 + t121 * t77, t11 * t28 + t19 * t13 - t2 * t47 + (-(-qJD(6) * t20 + t21) * t129 - t117 - (-qJD(6) * t6 + t16) * t50 + t5 * t105) * t77 + t121 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t35, 0, 0, 0, 0, 0, t82 - t91, t119 + t122; 0, 0, 0, 0, 0, -t81, -t130 * qJD(1), -t81 * t73, -t81 * t75, 0 (-t111 * qJD(4) - t51) * qJD(1), 0, 0, 0, 0, 0, -qJD(1) * t86 + t35, -qJD(1) * t45 - t36, 0, 0, 0, 0, 0, -t50 * t112 + t49 * t14 - t46 * t26 + (-t47 * t77 - t92 * t79) * t129, -t50 * t32 + t114 - t46 * t28 + (-t47 * t79 + t92 * t77) * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111 * t81, t90 * qJD(1) + t71, 0, 0, 0, 0, 0, -t54 + (t45 + t98) * qJD(5), -0.2e1 * t126, 0, 0, 0, 0, 0, t89 - t115, -t116 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t86, t45 ^ 2 - t86 ^ 2, 0, t54 + (t45 - t98) * qJD(5), 0, t49 * t101 - t40 * t45 (qJD(4) + t40) * t86, t28 * t94 + t118 (t13 - t131) * t79 + (-t129 * t28 - t14) * t77, -t116 - t128, t89 + t115, -t129 * t45, -pkin(5) * t14 - t1 * t45 - t123 * t79 - t8 * t26 + t83 * t77, -pkin(5) * t13 + t123 * t77 + t2 * t45 - t8 * t28 + t83 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t26, -t26 ^ 2 + t28 ^ 2, t13 + t131, t124 * t28 + t113, t34, t124 * t2 - t5 * t28 - t77 * t3 + t15, t1 * t124 - t77 * t16 + t5 * t26 - t79 * t3;];
tauc_reg  = t9;
