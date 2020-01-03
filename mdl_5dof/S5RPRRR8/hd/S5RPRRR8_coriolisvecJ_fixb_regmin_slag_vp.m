% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:06
% EndTime: 2019-12-31 19:06:09
% DurationCPUTime: 0.79s
% Computational Cost: add. (979->126), mult. (1587->194), div. (0->0), fcn. (944->6), ass. (0->105)
t61 = qJD(4) + qJD(5);
t70 = cos(qJ(3));
t141 = 0.2e1 * t70;
t65 = sin(qJ(5));
t66 = sin(qJ(4));
t68 = cos(qJ(5));
t69 = cos(qJ(4));
t38 = t65 * t69 + t68 * t66;
t135 = t38 * t61;
t124 = t135 * t61;
t107 = qJD(1) - qJD(3);
t10 = t135 * t107;
t140 = t107 ^ 2;
t129 = pkin(7) + pkin(8);
t67 = sin(qJ(3));
t72 = qJD(4) ^ 2;
t139 = (-t72 - t140) * t67;
t130 = qJD(5) - t61;
t116 = t68 * t69;
t138 = t61 * t116;
t108 = (qJD(1) * qJD(2));
t71 = -pkin(1) - pkin(2);
t50 = t71 * qJD(1) + qJD(2);
t137 = qJD(3) * t50 + t108;
t136 = qJD(4) * t107;
t109 = qJD(1) * qJ(2);
t94 = qJD(3) * t109;
t20 = t137 * t67 + t70 * t94;
t32 = t70 * t109 + t67 * t50;
t133 = -t107 * t32 - t20;
t85 = t70 * qJ(2) + t67 * t71;
t30 = t67 * qJD(2) + qJD(3) * t85;
t132 = t107 * t30 + t20;
t131 = -t67 * qJ(2) + t70 * t71;
t128 = t107 * pkin(3);
t127 = t69 * pkin(4);
t40 = -pkin(7) + t85;
t126 = pkin(8) - t40;
t120 = t65 * t66;
t17 = t61 * t120 - t138;
t125 = t17 * t61;
t105 = t107 * t120;
t26 = t107 * t116 - t105;
t28 = t38 * t107;
t123 = t28 * t26;
t98 = -t129 * t107 + t32;
t15 = t98 * t69;
t118 = t68 * t15;
t114 = t72 * t66;
t113 = t72 * t69;
t112 = t66 ^ 2 - t69 ^ 2;
t110 = qJD(4) * t66;
t106 = pkin(4) * t107 * t66;
t104 = pkin(4) * t110;
t103 = t107 * t110;
t100 = 2 * t108;
t58 = -pkin(3) - t127;
t14 = t98 * t66;
t12 = qJD(4) * pkin(4) - t14;
t99 = -pkin(4) * t61 - t12;
t97 = qJD(4) * t129;
t19 = t137 * t70 - t67 * t94;
t31 = -t67 * t109 + t70 * t50;
t22 = -t31 + t128;
t96 = t107 * t22 - t19;
t95 = qJD(4) * t126;
t92 = t69 * t103;
t39 = pkin(3) - t131;
t91 = -t32 + t104;
t9 = t61 * t105 - t138 * t107;
t90 = t28 * t17 + t9 * t38;
t13 = -pkin(4) * t103 + t20;
t16 = -t107 * t58 - t31;
t88 = -t13 * t38 + t16 * t17;
t37 = -t116 + t120;
t87 = -t13 * t37 - t135 * t16;
t86 = qJD(4) * t98;
t4 = t69 * t19 - t66 * t86;
t5 = -t66 * t19 - t69 * t86;
t84 = t16 * t28 - t65 * t4 + t68 * t5;
t83 = pkin(7) * t72 - t133;
t82 = t40 * t72 - t132;
t81 = qJD(4) * (t22 + t31 + t128);
t29 = t70 * qJD(2) + t131 * qJD(3);
t80 = qJD(4) * (-t107 * t39 - t22 - t29);
t79 = t136 * t141;
t77 = t37 * t61;
t76 = t16 * t26 + (t130 * t15 - t5) * t65;
t75 = -t38 * t10 - t135 * t28 - t17 * t26 + t9 * t37;
t73 = qJD(1) ^ 2;
t49 = t129 * t69;
t48 = t129 * t66;
t42 = t69 * t97;
t41 = t66 * t97;
t34 = t39 + t127;
t33 = t112 * t136;
t25 = t126 * t69;
t24 = t126 * t66;
t21 = t30 - t104;
t8 = -t66 * t29 + t69 * t95;
t7 = t69 * t29 + t66 * t95;
t6 = -t26 ^ 2 + t28 ^ 2;
t3 = -t28 * t61 + t10;
t2 = t26 * t61 + t9;
t1 = [0, 0, 0, 0, t100, qJ(2) * t100, 0, t132, t107 * t29 + t19, 0.2e1 * t92, -0.2e1 * t33, -t113, t114, 0, t66 * t80 - t69 * t82, t66 * t82 + t69 * t80, -t90, t75, t125, t124, 0, t21 * t26 - t34 * t10 + (-t65 * t7 + t68 * t8 + (-t24 * t65 + t25 * t68) * qJD(5)) * t61 + t87, -t21 * t28 + t34 * t9 - (t65 * t8 + t68 * t7 + (t24 * t68 + t25 * t65) * qJD(5)) * t61 + t88; 0, 0, 0, 0, -t73, -t73 * qJ(2), 0, -t67 * t140, -t70 * t140, 0, 0, 0, 0, 0, t69 * t139 + t66 * t79, -t66 * t139 + t69 * t79, 0, 0, 0, 0, 0, t10 * t141 + (t61 ^ 2 * t37 - t107 * t26) * t67, (-t107 * t77 - t9) * t70 + (t107 * t28 + t124) * t67; 0, 0, 0, 0, 0, 0, 0, t133, -t107 * t31 - t19, -0.2e1 * t92, 0.2e1 * t33, t113, -t114, 0, t66 * t81 - t69 * t83, t66 * t83 + t69 * t81, t90, -t75, -t125, -t124, 0, -t58 * t10 + (t65 * t41 - t68 * t42 + (t48 * t65 - t49 * t68) * qJD(5)) * t61 + t31 * t135 + t91 * t26 - t87, t58 * t9 - (-t68 * t41 - t65 * t42 + (-t48 * t68 - t49 * t65) * qJD(5)) * t61 - t31 * t77 - t91 * t28 - t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66 * t140 * t69, t112 * t140, 0, 0, 0, t96 * t66, t96 * t69, -t123, t6, t2, t3, 0, t26 * t106 - (t65 * t14 - t118) * t61 + (t99 * t65 - t118) * qJD(5) + t84, -t28 * t106 + (t99 * qJD(5) - t14 * t61 - t4) * t68 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, t6, t2, t3, 0, t84 + t130 * (-t65 * t12 - t118), (-t130 * t12 - t4) * t68 + t76;];
tauc_reg = t1;
