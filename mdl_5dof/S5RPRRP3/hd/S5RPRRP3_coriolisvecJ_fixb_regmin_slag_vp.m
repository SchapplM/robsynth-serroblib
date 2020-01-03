% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:47:42
% EndTime: 2020-01-03 11:47:45
% DurationCPUTime: 0.65s
% Computational Cost: add. (1002->127), mult. (2405->181), div. (0->0), fcn. (1560->6), ass. (0->93)
t118 = cos(qJ(4));
t70 = sin(qJ(4));
t71 = sin(qJ(3));
t72 = cos(qJ(3));
t50 = t118 * t71 + t70 * t72;
t103 = qJD(1) * t50;
t104 = t103 * qJ(5);
t60 = sin(pkin(8)) * pkin(1) + pkin(6);
t119 = pkin(7) + t60;
t87 = t119 * qJD(1);
t33 = t71 * qJD(2) + t72 * t87;
t28 = t70 * t33;
t106 = qJD(3) * pkin(3);
t32 = t72 * qJD(2) - t87 * t71;
t31 = t32 + t106;
t91 = t118 * t31 - t28;
t122 = t104 - t91;
t65 = qJD(3) + qJD(4);
t121 = t103 ^ 2;
t5 = t65 * pkin(4) - t122;
t120 = t5 + t122;
t95 = t118 * t72;
t86 = qJD(1) * t95;
t102 = qJD(1) * t71;
t96 = t70 * t102;
t42 = -t86 + t96;
t61 = -cos(pkin(8)) * pkin(1) - pkin(2);
t51 = -t72 * pkin(3) + t61;
t46 = t51 * qJD(1);
t20 = t42 * pkin(4) + qJD(5) + t46;
t117 = t20 * t103;
t113 = t70 * t71;
t85 = t65 * t113;
t88 = t118 * qJD(4);
t23 = -qJD(3) * t95 - t72 * t88 + t85;
t116 = t23 * t65;
t115 = t103 * t42;
t114 = t46 * t103;
t73 = qJD(3) ^ 2;
t112 = t73 * t71;
t111 = t73 * t72;
t24 = t65 * t50;
t18 = t24 * qJD(1);
t110 = -t50 * t18 + t23 * t42;
t109 = t118 * t32 - t28;
t108 = t65 * t86;
t107 = t71 ^ 2 - t72 ^ 2;
t105 = t42 * qJ(5);
t53 = qJD(1) * t61;
t101 = qJD(4) * t70;
t100 = t53 * qJD(1);
t98 = qJD(1) * qJD(3);
t97 = t71 * t106;
t30 = t118 * t33;
t93 = t71 * t98;
t59 = pkin(3) * t93;
t94 = t18 * pkin(4) + t59;
t26 = t32 * qJD(3);
t27 = t33 * qJD(3);
t92 = -t118 * t27 - t70 * t26;
t90 = -t70 * t32 - t30;
t89 = qJD(3) * t119;
t17 = qJD(1) * t85 - t108;
t49 = -t95 + t113;
t84 = t103 * t24 - t49 * t17;
t82 = 0.2e1 * qJD(3) * t53;
t81 = -t70 * t31 - t30;
t47 = t119 * t71;
t48 = t119 * t72;
t80 = -t118 * t48 + t70 * t47;
t39 = t71 * t89;
t40 = t72 * t89;
t79 = -t101 * t48 - t118 * t39 - t70 * t40 - t47 * t88;
t78 = qJD(4) * t81 + t92;
t77 = qJD(4) * t80 - t118 * t40 + t70 * t39;
t76 = -t33 * t101 + t118 * t26 - t70 * t27 + t31 * t88;
t75 = t46 * t42 - t76;
t74 = qJD(1) ^ 2;
t63 = pkin(3) * t118 + pkin(4);
t41 = t42 ^ 2;
t19 = t24 * t65;
t15 = -t41 + t121;
t13 = -t49 * qJ(5) - t80;
t12 = -t50 * qJ(5) - t118 * t47 - t70 * t48;
t10 = t108 + (t42 - t96) * t65;
t9 = -t104 + t109;
t8 = t90 + t105;
t7 = -t81 - t105;
t4 = t23 * qJ(5) - t50 * qJD(5) + t77;
t3 = -t24 * qJ(5) - t49 * qJD(5) + t79;
t2 = t17 * qJ(5) - qJD(5) * t103 + t78;
t1 = -t18 * qJ(5) - t42 * qJD(5) + t76;
t6 = [0, 0, 0, 0, 0.2e1 * t72 * t93, -0.2e1 * t107 * t98, t111, -t112, 0, -t111 * t60 + t71 * t82, t112 * t60 + t72 * t82, -t103 * t23 - t17 * t50, -t84 + t110, -t116, -t19, 0, t51 * t18 + t46 * t24 + t42 * t97 + t59 * t49 + t65 * t77, 0.2e1 * t103 * t97 - t51 * t17 - t46 * t23 - t79 * t65, -t1 * t49 - t103 * t4 + t12 * t17 - t13 * t18 - t2 * t50 + t5 * t23 - t7 * t24 - t3 * t42, t1 * t13 + t7 * t3 + t2 * t12 + t5 * t4 + t94 * (t49 * pkin(4) + t51) + t20 * (t24 * pkin(4) + t97); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t111, 0, 0, 0, 0, 0, -t19, t116, t84 + t110, t1 * t50 - t2 * t49 - t7 * t23 - t5 * t24; 0, 0, 0, 0, -t71 * t74 * t72, t107 * t74, 0, 0, 0, -t71 * t100, -t72 * t100, t115, t15, t10, 0, 0, -pkin(3) * t42 * t102 - t114 - t90 * t65 + (-t30 + (-pkin(3) * t65 - t31) * t70) * qJD(4) + t92, t109 * t65 + (-t102 * t103 - t65 * t88) * pkin(3) + t75, t63 * t17 + (t7 + t8) * t103 + (-t5 + t9) * t42 + (-t18 * t70 + (t103 * t70 - t118 * t42) * qJD(4)) * pkin(3), -pkin(4) * t117 + t2 * t63 - t5 * t8 - t7 * t9 + (-t20 * t102 + t1 * t70 + (t118 * t7 - t5 * t70) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t15, t10, 0, 0, -t65 * t81 - t114 + t78, t65 * t91 + t75, pkin(4) * t17 - t120 * t42, t120 * t7 + (t2 - t117) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41 - t121, t103 * t5 + t7 * t42 + t94;];
tauc_reg = t6;
