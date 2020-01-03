% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:06
% EndTime: 2019-12-31 17:17:08
% DurationCPUTime: 0.68s
% Computational Cost: add. (1151->163), mult. (3115->210), div. (0->0), fcn. (1986->4), ass. (0->103)
t106 = qJD(2) * pkin(2);
t71 = sin(qJ(2));
t99 = t71 * t106;
t129 = 0.2e1 * t99;
t67 = qJD(2) + qJD(3);
t101 = (qJD(1) * qJD(2));
t128 = -2 * t101;
t70 = sin(qJ(3));
t72 = cos(qJ(2));
t113 = t70 * t72;
t123 = cos(qJ(3));
t48 = t123 * t71 + t113;
t105 = qJD(1) * t48;
t127 = t105 ^ 2;
t126 = -pkin(6) - pkin(5);
t54 = t126 * t71;
t50 = qJD(1) * t54;
t44 = t50 + t106;
t55 = t126 * t72;
t52 = qJD(1) * t55;
t96 = t123 * t52;
t26 = t70 * t44 - t96;
t94 = qJD(2) * t126;
t84 = qJD(1) * t94;
t45 = t71 * t84;
t46 = t72 * t84;
t92 = -t123 * t46 + t70 * t45;
t7 = t26 * qJD(3) + t92;
t79 = t123 * t54 + t70 * t55;
t125 = t7 * t79;
t51 = t71 * t94;
t9 = t79 * qJD(3) + t94 * t113 + t123 * t51;
t124 = t9 * t67;
t32 = -t123 * t55 + t70 * t54;
t95 = t123 * t72;
t87 = qJD(2) * t95;
t10 = t32 * qJD(3) - t126 * t87 + t70 * t51;
t122 = t10 * t67;
t86 = qJD(1) * t95;
t104 = qJD(1) * t71;
t97 = t70 * t104;
t38 = -t86 + t97;
t65 = -t72 * pkin(2) - pkin(1);
t53 = qJD(1) * t65;
t16 = t38 * pkin(3) - qJ(4) * t105 + t53;
t121 = t16 * t105;
t115 = t70 * t52;
t25 = t123 * t44 + t115;
t120 = t25 * t67;
t119 = t105 * t38;
t117 = t53 * t105;
t30 = t67 * t48;
t116 = t67 * t30;
t114 = t70 * t71;
t74 = qJD(1) ^ 2;
t112 = t72 * t74;
t73 = qJD(2) ^ 2;
t111 = t73 * t71;
t110 = t73 * t72;
t28 = t123 * t50 + t115;
t91 = t123 * qJD(3);
t109 = pkin(2) * t91 + qJD(4) - t28;
t108 = t67 * t86;
t107 = t71 ^ 2 - t72 ^ 2;
t103 = qJD(3) * t70;
t102 = qJD(4) - t25;
t100 = t71 * t112;
t98 = pkin(2) * t104;
t14 = t38 ^ 2 - t127;
t93 = t71 * t101;
t90 = -t52 * t103 - t123 * t45 - t44 * t91 - t70 * t46;
t89 = pkin(1) * t128;
t88 = t72 * t93;
t27 = t70 * t50 - t96;
t85 = pkin(2) * t103 - t27;
t83 = t67 * t114;
t22 = pkin(3) * t105 + t38 * qJ(4);
t21 = t30 * qJD(1);
t47 = -t95 + t114;
t82 = t21 * t47 + t38 * t30;
t81 = -t16 * t38 - t90;
t80 = t53 * t38 + t90;
t20 = qJD(1) * t83 - t108;
t78 = t10 * t105 + t20 * t79 - t32 * t21 - t9 * t38 + t7 * t48;
t29 = -t72 * t91 + t83 - t87;
t77 = -t105 * t30 + t20 * t47 - t48 * t21 + t29 * t38;
t76 = -t92 + (-qJD(3) + t67) * t26;
t75 = t27 * t67 + (t96 + (-pkin(2) * t67 - t44) * t70) * qJD(3) - t92;
t66 = t67 * qJD(4);
t64 = -t123 * pkin(2) - pkin(3);
t62 = t70 * pkin(2) + qJ(4);
t24 = t29 * t67;
t23 = t47 * pkin(3) - t48 * qJ(4) + t65;
t19 = t67 * qJ(4) + t26;
t18 = t22 + t98;
t17 = -t67 * pkin(3) + t102;
t12 = -t67 * t105 + t21;
t11 = t108 + (t38 - t97) * t67;
t4 = t66 - t90;
t3 = t30 * pkin(3) + t29 * qJ(4) - t48 * qJD(4) + t99;
t2 = -t105 * t29 - t20 * t48;
t1 = pkin(2) * t93 + t21 * pkin(3) + t20 * qJ(4) - qJD(4) * t105;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t88, t107 * t128, t110, -0.2e1 * t88, -t111, 0, -pkin(5) * t110 + t71 * t89, pkin(5) * t111 + t72 * t89, 0, 0, t2, t77, -t24, t82, -t116, 0, -t122 + t65 * t21 + t53 * t30 + (qJD(1) * t47 + t38) * t99, t105 * t129 - t65 * t20 - t53 * t29 - t124, t25 * t29 - t26 * t30 + t47 * t90 + t78, -t25 * t10 + t53 * t129 + t26 * t9 - t32 * t90 - t125, t2, -t24, -t77, 0, t116, t82, t1 * t47 + t16 * t30 + t23 * t21 + t3 * t38 - t122, -t17 * t29 - t19 * t30 - t4 * t47 + t78, -t1 * t48 - t105 * t3 + t16 * t29 + t23 * t20 + t124, t1 * t23 + t17 * t10 + t16 * t3 + t19 * t9 + t4 * t32 - t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t107 * t74, 0, t100, 0, 0, t74 * pkin(1) * t71, pkin(1) * t112, 0, 0, t119, -t14, t11, -t119, 0, 0, -t38 * t98 - t117 + t75, t28 * t67 + (-t104 * t105 - t67 * t91) * pkin(2) + t80, (t26 - t27) * t105 + (-t25 + t28) * t38 + (t123 * t20 - t21 * t70 + (t105 * t70 - t123 * t38) * qJD(3)) * pkin(2), t25 * t27 - t26 * t28 + (-t53 * t104 - t123 * t7 - t90 * t70 + (t123 * t26 - t25 * t70) * qJD(3)) * pkin(2), t119, t11, t14, 0, t12, -t119, -t18 * t38 - t121 + t75, -t64 * t20 - t62 * t21 + (t19 + t85) * t105 + (t17 - t109) * t38, t105 * t18 + t109 * t67 + t66 + t81, t109 * t19 - t16 * t18 + t85 * t17 + t4 * t62 + t7 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t14, t11, -t119, 0, 0, t76 - t117, t80 + t120, 0, 0, t119, t11, t14, 0, t12, -t119, -t22 * t38 - t121 + t76, pkin(3) * t20 - t21 * qJ(4) + (t19 - t26) * t105 + (t17 - t102) * t38, t105 * t22 - t120 + 0.2e1 * t66 + t81, -t7 * pkin(3) + t4 * qJ(4) + t102 * t19 - t16 * t22 - t17 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t11, -t67 ^ 2 - t127, -t19 * t67 + t121 + t7;];
tauc_reg = t5;
