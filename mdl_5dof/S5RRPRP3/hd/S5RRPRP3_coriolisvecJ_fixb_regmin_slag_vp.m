% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:13
% EndTime: 2019-12-31 19:51:16
% DurationCPUTime: 0.70s
% Computational Cost: add. (1479->161), mult. (2566->203), div. (0->0), fcn. (1662->6), ass. (0->107)
t86 = sin(pkin(8));
t87 = cos(pkin(8));
t126 = t86 ^ 2 + t87 ^ 2;
t124 = pkin(1) * qJD(2);
t112 = qJD(1) * t124;
t85 = qJD(1) + qJD(2);
t91 = cos(qJ(2));
t64 = t85 * qJD(3) + t91 * t112;
t145 = t126 * t64;
t90 = cos(qJ(4));
t132 = t90 * t87;
t88 = sin(qJ(4));
t135 = t88 * t86;
t65 = -t132 + t135;
t66 = t90 * t86 + t88 * t87;
t125 = pkin(1) * qJD(1);
t98 = -t91 * t125 + qJD(3);
t49 = t66 * t85;
t89 = sin(qJ(2));
t144 = t87 * t89;
t71 = (-pkin(7) - qJ(3)) * t86;
t81 = t87 * pkin(7);
t72 = t87 * qJ(3) + t81;
t96 = t90 * t71 - t88 * t72;
t131 = t96 * qJD(4) - t98 * t65;
t143 = t131 * qJD(4);
t142 = t49 ^ 2;
t123 = qJD(4) * t88;
t114 = t86 * t123;
t122 = qJD(4) * t90;
t113 = t87 * t122;
t69 = t85 * t113;
t43 = t85 * t114 - t69;
t61 = t66 * qJD(4);
t44 = t85 * t61;
t77 = t89 * t112;
t94 = t44 * pkin(4) + t43 * qJ(5) + t77;
t10 = -t49 * qJD(5) + t94;
t80 = -t87 * pkin(3) - pkin(2);
t46 = t80 * t85 + t98;
t118 = t85 * t135;
t47 = -t85 * t132 + t118;
t16 = t47 * pkin(4) - t49 * qJ(5) + t46;
t141 = t10 * t65 + t16 * t61;
t60 = -t113 + t114;
t140 = -t10 * t66 + t16 * t60;
t139 = t91 * pkin(1);
t138 = t16 * t49;
t137 = t49 * t47;
t116 = t89 * t125;
t68 = t85 * qJ(3) + t116;
t111 = pkin(7) * t85 + t68;
t40 = t111 * t87;
t136 = t88 * t40;
t38 = t88 * t71 + t90 * t72;
t130 = t38 * qJD(4) + t98 * t66;
t129 = t46 * t61 + t65 * t77;
t128 = -t46 * t60 + t66 * t77;
t75 = t91 * t124 + qJD(3);
t79 = t89 * pkin(1) + qJ(3);
t62 = (-pkin(7) - t79) * t86;
t63 = t87 * t79 + t81;
t97 = t90 * t62 - t88 * t63;
t14 = t97 * qJD(4) - t65 * t75;
t121 = t14 * qJD(4);
t32 = t88 * t62 + t90 * t63;
t15 = t32 * qJD(4) + t66 * t75;
t120 = t15 * qJD(4);
t39 = t111 * t86;
t20 = -t90 * t39 - t136;
t119 = qJD(5) - t20;
t117 = t89 * t124;
t110 = t126 * t91;
t108 = t126 * t75;
t17 = -qJD(4) * pkin(4) + t119;
t21 = -t88 * t39 + t90 * t40;
t18 = qJD(4) * qJ(5) + t21;
t95 = t39 * t122 + t65 * t64;
t4 = (qJD(5) - t136) * qJD(4) - t95;
t5 = t40 * t122 - t39 * t123 + t66 * t64;
t107 = -t17 * t60 - t18 * t61 - t4 * t65 + t5 * t66;
t106 = t20 + t136;
t105 = t130 * qJD(4);
t104 = t126 * qJD(3);
t103 = t85 * t117;
t102 = t85 * t116;
t22 = t61 * pkin(4) + t60 * qJ(5) - t66 * qJD(5);
t101 = -t22 + t116;
t100 = (-qJD(2) + t85) * t125;
t99 = (-qJD(1) - t85) * t124;
t93 = t21 * qJD(4) - t5;
t33 = t65 * pkin(4) - t66 * qJ(5) + t80;
t92 = 0.2e1 * t49 * qJD(4);
t73 = t86 * t77;
t70 = t80 - t139;
t67 = -t85 * pkin(2) + t98;
t57 = t61 * qJD(4);
t56 = t60 * qJD(4);
t45 = t47 ^ 2;
t28 = t33 - t139;
t27 = t49 * pkin(4) + t47 * qJ(5);
t26 = t69 + (t47 - t118) * qJD(4);
t25 = -t69 + (t47 + t118) * qJD(4);
t19 = t22 + t117;
t13 = -t43 * t66 - t49 * t60;
t1 = t43 * t65 - t66 * t44 + t60 * t47 - t49 * t61;
t2 = [0, 0, 0, 0, -t77 - t103, t91 * t99, t99 * t144, t103 * t86 + t73, t85 * t108 + t145, t68 * t108 + t79 * t145 + (t67 + (-pkin(2) - t139) * qJD(1)) * t117, t13, t1, -t56, -t57, 0, t47 * t117 + t70 * t44 - t120 + t129, t49 * t117 - t70 * t43 - t121 + t128, t19 * t47 + t28 * t44 - t120 + t141, -t14 * t47 + t15 * t49 - t32 * t44 + t43 * t97 + t107, -t19 * t49 + t28 * t43 + t121 + t140, t10 * t28 + t18 * t14 + t17 * t15 + t16 * t19 + t4 * t32 - t5 * t97; 0, 0, 0, 0, -t77 + t102, t91 * t100, t100 * t144, -t102 * t86 + t73, (-t110 * t125 + t104) * t85 + t145, t68 * t104 + qJ(3) * t145 + ((-pkin(2) * qJD(2) - t67) * t89 - t68 * t110) * t125, t13, t1, -t56, -t57, 0, -t47 * t116 + t80 * t44 - t105 + t129, -t49 * t116 - t80 * t43 + t128 - t143, -t101 * t47 + t33 * t44 - t105 + t141, t130 * t49 - t131 * t47 - t38 * t44 + t43 * t96 + t107, t101 * t49 + t33 * t43 + t140 + t143, t10 * t33 - t101 * t16 + t130 * t17 + t131 * t18 + t4 * t38 - t5 * t96; 0, 0, 0, 0, 0, 0, 0, 0, -t126 * t85 ^ 2, -t126 * t85 * t68 + t77, 0, 0, 0, 0, 0, t92, -t25, t92, -t45 - t142, t25, t18 * t47 + (-qJD(5) - t17) * t49 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, -t45 + t142, t26, 0, 0, -t46 * t49 + t93, t106 * qJD(4) + t46 * t47 + t95, -t27 * t47 - t138 + t93, pkin(4) * t43 - t44 * qJ(5) + (t18 - t21) * t49 + (t17 - t119) * t47, -t16 * t47 + t27 * t49 + (0.2e1 * qJD(5) - t106) * qJD(4) - t95, -t5 * pkin(4) + t4 * qJ(5) + t119 * t18 - t16 * t27 - t17 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t26, -qJD(4) ^ 2 - t142, -t18 * qJD(4) + t138 + t5;];
tauc_reg = t2;
