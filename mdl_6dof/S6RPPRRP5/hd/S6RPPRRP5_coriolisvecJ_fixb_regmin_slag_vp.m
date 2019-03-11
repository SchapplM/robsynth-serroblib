% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:47
% EndTime: 2019-03-09 02:08:50
% DurationCPUTime: 1.16s
% Computational Cost: add. (1365->224), mult. (2697->319), div. (0->0), fcn. (1466->4), ass. (0->123)
t71 = cos(qJ(4));
t120 = qJD(1) * t71;
t70 = cos(qJ(5));
t100 = t70 * t120;
t68 = sin(qJ(5));
t44 = t68 * qJD(4) + t100;
t69 = sin(qJ(4));
t112 = t69 * qJD(1);
t58 = qJD(5) + t112;
t137 = t44 * t58;
t108 = qJD(1) * qJD(4);
t97 = t69 * t108;
t17 = t44 * qJD(5) - t68 * t97;
t151 = t17 - t137;
t117 = qJD(4) * t71;
t66 = -pkin(7) + qJ(2);
t150 = qJD(2) * t69 + t66 * t117;
t111 = t70 * qJD(4);
t42 = t68 * t120 - t111;
t138 = t42 * t58;
t114 = qJD(5) * t71;
t102 = t68 * t114;
t77 = t69 * t111 + t102;
t16 = t77 * qJD(1) - qJD(5) * t111;
t149 = -t16 + t138;
t67 = pkin(1) + qJ(3);
t148 = qJD(1) * t67;
t147 = t44 ^ 2;
t48 = t69 * pkin(4) - t71 * pkin(8) + t67;
t25 = t48 * qJD(1) - qJD(2);
t60 = qJD(1) * qJ(2) + qJD(3);
t55 = -pkin(7) * qJD(1) + t60;
t47 = t69 * t55;
t32 = qJD(4) * pkin(8) + t47;
t10 = t70 * t25 - t68 * t32;
t6 = -t44 * qJ(6) + t10;
t5 = t58 * pkin(5) + t6;
t146 = t5 - t6;
t121 = t70 * qJ(6);
t129 = t71 * t55;
t87 = pkin(4) * t71 + pkin(8) * t69;
t46 = t87 * qJD(1);
t88 = -t68 * t129 + t70 * t46;
t128 = -qJ(6) - pkin(8);
t94 = qJD(5) * t128;
t145 = -t68 * qJD(6) + t70 * t94 - (pkin(5) * t71 + t69 * t121) * qJD(1) - t88;
t144 = t16 * t68;
t143 = t16 * t69;
t142 = t17 * t69;
t109 = qJD(1) * qJD(2);
t118 = qJD(4) * t69;
t23 = -t71 * t109 + t55 * t118;
t141 = t23 * t68;
t140 = t23 * t70;
t33 = -qJD(4) * pkin(4) - t129;
t139 = t33 * t70;
t136 = t44 * t71;
t135 = t58 * t68;
t134 = t58 * t70;
t133 = t66 * t68;
t132 = t68 * t69;
t131 = t69 * t70;
t130 = t71 * t16;
t126 = t70 * t129 + t68 * t46;
t127 = t70 * qJD(6) - t126 + (-qJ(6) * t112 + t94) * t68;
t125 = t66 * t131 + t68 * t48;
t65 = t71 ^ 2;
t124 = t69 ^ 2 - t65;
t72 = qJD(4) ^ 2;
t73 = qJD(1) ^ 2;
t123 = -t72 - t73;
t122 = qJ(6) * t71;
t116 = qJD(5) * t68;
t115 = qJD(5) * t70;
t113 = t33 * qJD(5);
t56 = -qJD(2) + t148;
t110 = qJD(2) - t56;
t107 = t58 * t131;
t106 = t58 * t132;
t24 = t69 * t109 + t55 * t117;
t40 = t87 * qJD(4) + qJD(3);
t26 = t40 * qJD(1);
t105 = -t25 * t115 - t70 * t24 - t68 * t26;
t103 = t58 * t116;
t101 = t70 * t114;
t62 = 0.2e1 * t109;
t99 = 0.2e1 * qJD(3) * qJD(1);
t98 = pkin(5) * t68 - t66;
t96 = t71 * t108;
t95 = t58 * t66 + t32;
t93 = t48 * t115 + t150 * t70 + t68 * t40;
t92 = t58 + t112;
t91 = t110 * qJD(1);
t90 = -0.2e1 * t96;
t89 = qJD(5) * t69 + qJD(1);
t86 = t58 * t115 + t68 * t96;
t11 = t68 * t25 + t70 * t32;
t7 = -t42 * qJ(6) + t11;
t85 = t5 * t70 + t68 * t7;
t84 = t5 * t68 - t7 * t70;
t83 = qJD(2) + t56 + t148;
t82 = qJD(1) * t65 - t58 * t69;
t80 = -t66 * t72 + t99;
t79 = -pkin(8) * t117 + t33 * t69;
t78 = -t32 * t116 - t105;
t8 = t17 * pkin(5) + t23;
t76 = -qJD(6) * t71 + (qJ(6) * qJD(4) - qJD(5) * t66) * t69;
t22 = t70 * t26;
t75 = -t11 * qJD(5) - t68 * t24 + t22;
t1 = pkin(5) * t96 + t16 * qJ(6) - t44 * qJD(6) + t75;
t2 = -t17 * qJ(6) - t42 * qJD(6) + t78;
t74 = t84 * qJD(5) - t1 * t70 - t2 * t68;
t51 = t128 * t70;
t50 = t128 * t68;
t39 = t42 ^ 2;
t37 = t70 * t48;
t30 = t70 * t40;
t15 = -t68 * t122 + t125;
t14 = t42 * pkin(5) + qJD(6) + t33;
t13 = -t71 * t121 + t37 + (pkin(5) - t133) * t69;
t4 = -qJ(6) * t101 + t76 * t68 + t93;
t3 = pkin(5) * t117 + t30 + t76 * t70 + ((-t48 + t122) * qJD(5) - t150) * t68;
t9 = [0, 0, 0, 0, t62, qJ(2) * t62, t62, t99, t60 * qJD(2) + t56 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t67) * qJD(1), t69 * t90, 0.2e1 * t124 * t108, -t72 * t69, -t72 * t71, 0, t83 * t117 + t80 * t69, -t83 * t118 + t71 * t80, -t70 * t130 - t44 * t77 (t42 * t70 + t44 * t68) * t118 + (t144 - t17 * t70 + (t42 * t68 - t44 * t70) * qJD(5)) * t71, -t58 * t102 - t143 + (t70 * t82 + t136) * qJD(4), -t58 * t101 - t142 + (-t71 * t42 - t68 * t82) * qJD(4), t92 * t117 (-t48 * t116 + t30) * t58 + (qJD(4) * t66 * t42 + t22 - t95 * t115 + (-qJD(2) * t58 - t33 * qJD(4) - qJD(5) * t25 - t24) * t68) * t69 + (t70 * t113 - qJD(2) * t42 - t66 * t17 + t141 + (-t58 * t133 + (-t66 * t132 + t37) * qJD(1) + t10) * qJD(4)) * t71, -t93 * t58 + (t95 * t116 + (t44 * t66 - t139) * qJD(4) + t105) * t69 + (-t68 * t113 - qJD(2) * t44 + t66 * t16 + t140 + (-t125 * qJD(1) - t11) * qJD(4)) * t71, t85 * t118 + t13 * t16 - t15 * t17 - t3 * t44 - t4 * t42 + t71 * t74, t1 * t13 + t2 * t15 + t5 * t3 + t7 * t4 - t14 * t98 * t118 + (t8 * t98 + t14 * (pkin(5) * t115 - qJD(2))) * t71; 0, 0, 0, 0, -t73, -t73 * qJ(2), -t73, 0 (-qJD(3) - t60) * qJD(1), 0, 0, 0, 0, 0, t90, 0.2e1 * t97, 0, 0, 0, 0, 0, t103 + (t106 + (t42 - t111) * t71) * qJD(1) (t107 + t136) * qJD(1) + t86, t149 * t70 + t151 * t68 (t14 * t71 + t69 * t84) * qJD(1) + t74; 0, 0, 0, 0, 0, 0, 0, -t73, t91, 0, 0, 0, 0, 0, t123 * t69, t123 * t71, 0, 0, 0, 0, 0, -t71 * t17 - t89 * t134 + (-t68 * t71 * t92 + t69 * t42) * qJD(4), t130 + t89 * t135 + (-t71 * t134 + (t44 - t100) * t69) * qJD(4) (-t42 * t117 + t44 * t89 - t142) * t70 + (t44 * t117 + t42 * t89 - t143) * t68, -t85 * qJD(1) + (-qJD(4) * t84 - t8) * t71 + (qJD(4) * t14 - t85 * qJD(5) - t1 * t68 + t2 * t70) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t73 * t69, -t124 * t73, 0, 0, 0, t71 * t91, -t110 * t112, t44 * t134 - t144 (-t16 - t138) * t70 + (-t17 - t137) * t68 (t107 - t136) * qJD(1) + t86, -t103 + (-t106 + (t42 + t111) * t71) * qJD(1), -t58 * t120, -pkin(4) * t17 - t140 - t88 * t58 - t42 * t47 + (-pkin(8) * t134 + t33 * t68) * qJD(5) + (-t10 * t71 + t68 * t79) * qJD(1), pkin(4) * t16 + t141 + t126 * t58 - t44 * t47 + (pkin(8) * t135 + t139) * qJD(5) + (t11 * t71 + t70 * t79) * qJD(1), t50 * t16 + t51 * t17 - t145 * t44 - t127 * t42 + (-t5 * t58 + t2) * t70 + (-t58 * t7 - t1) * t68, -t2 * t51 + t1 * t50 + t8 * (-t70 * pkin(5) - pkin(4)) + t127 * t7 + t145 * t5 + (pkin(5) * t135 - t47) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t42, -t39 + t147, t149, -t151, t96, t11 * t58 - t33 * t44 + t75, t10 * t58 + t33 * t42 - t78, pkin(5) * t16 - t146 * t42, t146 * t7 + (-t14 * t44 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39 - t147, t7 * t42 + t5 * t44 + t8;];
tauc_reg  = t9;
