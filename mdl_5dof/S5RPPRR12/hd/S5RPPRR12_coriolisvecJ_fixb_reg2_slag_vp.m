% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR12_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:26
% DurationCPUTime: 1.44s
% Computational Cost: add. (2525->219), mult. (5495->297), div. (0->0), fcn. (3757->6), ass. (0->119)
t137 = sin(qJ(4));
t73 = sin(pkin(8));
t105 = t137 * t73;
t74 = cos(pkin(8));
t78 = cos(qJ(4));
t120 = t74 * t78;
t52 = -t105 + t120;
t51 = t137 * t74 + t78 * t73;
t82 = qJD(1) * t51;
t145 = qJD(5) + t82;
t83 = t51 * qJD(3);
t75 = -pkin(1) - qJ(3);
t142 = t75 * qJD(1);
t58 = qJD(2) + t142;
t104 = -pkin(6) * qJD(1) + t58;
t39 = t104 * t73;
t40 = t104 * t74;
t85 = t137 * t39 - t78 * t40;
t11 = -qJD(1) * t83 - t85 * qJD(4);
t48 = t51 * qJD(4);
t37 = qJD(1) * t48;
t108 = qJD(4) * t120;
t99 = qJD(1) * t105;
t57 = qJD(4) * t99;
t38 = qJD(1) * t108 - t57;
t71 = qJD(1) * qJD(2);
t22 = pkin(4) * t38 + pkin(7) * t37 + t71;
t24 = t137 * t40 + t78 * t39;
t20 = qJD(4) * pkin(7) + t24;
t109 = qJD(1) * t120;
t46 = -t99 + t109;
t72 = qJD(1) * qJ(2);
t66 = qJD(3) + t72;
t68 = t73 * pkin(3);
t55 = qJD(1) * t68 + t66;
t21 = pkin(4) * t82 - pkin(7) * t46 + t55;
t76 = sin(qJ(5));
t77 = cos(qJ(5));
t6 = t20 * t77 + t21 * t76;
t2 = -qJD(5) * t6 - t76 * t11 + t77 * t22;
t140 = t145 * t6 + t2;
t93 = t20 * t76 - t21 * t77;
t1 = -t93 * qJD(5) + t77 * t11 + t76 * t22;
t95 = t145 * t93 + t1;
t146 = t52 * qJD(3);
t101 = t77 * t145;
t118 = t76 * t38;
t144 = -t101 * t145 - t118;
t34 = qJD(4) * t76 + t46 * t77;
t115 = qJD(5) * t34;
t15 = -t76 * t37 + t115;
t116 = t73 ^ 2 + t74 ^ 2;
t141 = t116 * qJD(3);
t139 = t46 ^ 2;
t106 = 0.2e1 * t71;
t138 = -pkin(6) + t75;
t88 = t146 * qJD(1);
t12 = t24 * qJD(4) + t88;
t53 = t138 * t73;
t54 = t138 * t74;
t29 = t137 * t53 - t78 * t54;
t136 = t12 * t29;
t135 = t12 * t52;
t134 = t12 * t76;
t111 = t77 * qJD(4);
t114 = qJD(5) * t76;
t14 = -qJD(5) * t111 + t46 * t114 + t77 * t37;
t133 = t14 * t76;
t132 = t15 * t77;
t32 = t46 * t76 - t111;
t131 = t32 * t82;
t130 = t32 * t46;
t129 = t32 * t76;
t128 = t32 * t77;
t127 = t34 * t32;
t126 = t34 * t46;
t125 = t34 * t76;
t124 = t34 * t77;
t123 = t38 * t51;
t122 = t46 * t82;
t121 = t52 * t77;
t13 = t76 * t15;
t36 = t77 * t38;
t113 = qJD(5) * t77;
t117 = -t32 * t113 - t13;
t62 = qJ(2) + t68;
t49 = -qJD(4) * t105 + t108;
t112 = t49 * qJD(4);
t103 = qJD(1) * t116;
t102 = t76 * t145;
t100 = qJD(5) * t51 + qJD(1);
t97 = t6 * t76 - t77 * t93;
t96 = -t6 * t77 - t76 * t93;
t27 = pkin(4) * t51 - pkin(7) * t52 + t62;
t30 = t137 * t54 + t78 * t53;
t9 = t27 * t77 - t30 * t76;
t10 = t27 * t76 + t30 * t77;
t92 = t124 + t129;
t91 = -t37 * t52 - t46 * t48;
t90 = t49 * t82 + t123;
t89 = t36 + (-t76 * t82 - t114) * t145;
t87 = t52 * t113 - t48 * t76;
t86 = -t52 * t114 - t48 * t77;
t19 = -qJD(4) * pkin(4) + t85;
t84 = -pkin(7) * t38 + t145 * t19;
t81 = t11 * t51 + t24 * t49 + t48 * t85 - t135;
t80 = -t97 * qJD(5) + t1 * t77 - t2 * t76;
t79 = qJD(1) ^ 2;
t43 = t82 ^ 2;
t41 = t48 * qJD(4);
t28 = pkin(4) * t46 + pkin(7) * t82;
t25 = pkin(4) * t49 + pkin(7) * t48 + qJD(2);
t17 = t30 * qJD(4) + t146;
t16 = -t29 * qJD(4) - t83;
t8 = t28 * t76 - t77 * t85;
t7 = t28 * t77 + t76 * t85;
t4 = -t10 * qJD(5) - t76 * t16 + t77 * t25;
t3 = t9 * qJD(5) + t77 * t16 + t76 * t25;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, qJ(2) * t106, 0, 0, 0, 0, 0, 0, t73 * t106, t74 * t106, 0.2e1 * qJD(3) * t103, (t66 + t72) * qJD(2) + (-t58 - t142) * t141, t91, t37 * t51 - t38 * t52 - t46 * t49 + t48 * t82, -t41, t90, -t112, 0, 0.2e1 * t82 * qJD(2) - t17 * qJD(4) + t62 * t38 + t55 * t49, -t16 * qJD(4) - t62 * t37 - t55 * t48 + (qJD(1) * t52 + t46) * qJD(2), -t16 * t82 + t17 * t46 - t29 * t37 - t30 * t38 - t81, t11 * t30 + t136 + t24 * t16 + t85 * t17 + (qJD(1) * t62 + t55) * qJD(2), -t14 * t121 + t34 * t86, (t125 + t128) * t48 + (t133 - t132 + (-t124 + t129) * qJD(5)) * t52, -t14 * t51 + t145 * t86 + t34 * t49 + t52 * t36, t52 * t13 + t32 * t87, -t52 * t118 - t145 * t87 - t15 * t51 - t32 * t49, t145 * t49 + t123, t52 * t134 + t145 * t4 + t29 * t15 + t17 * t32 + t19 * t87 + t2 * t51 + t9 * t38 - t49 * t93, -t1 * t51 - t10 * t38 + t12 * t121 - t29 * t14 - t145 * t3 + t17 * t34 + t19 * t86 - t6 * t49, -t10 * t15 + t14 * t9 - t3 * t32 - t34 * t4 + t97 * t48 + (qJD(5) * t96 - t1 * t76 - t2 * t77) * t52, t1 * t10 + t17 * t19 + t2 * t9 + t3 * t6 - t4 * t93 + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t79 * qJ(2), 0, 0, 0, 0, 0, 0, -t79 * t73, -t79 * t74, 0, (-t66 - t141) * qJD(1), 0, 0, 0, 0, 0, 0, -qJD(1) * t82 - t41, -qJD(1) * t46 - t112, -t90 - t91, -qJD(1) * t55 + t81, 0, 0, 0, 0, 0, 0, -t51 * t118 - t52 * t15 + t48 * t32 + (-t100 * t77 - t49 * t76) * t145, -t51 * t36 + t52 * t14 + t48 * t34 + (t100 * t76 - t49 * t77) * t145, (t125 - t128) * t49 + t92 * qJD(1) + (qJD(5) * t92 - t132 - t133) * t51, -qJD(1) * t97 + t19 * t48 - t49 * t96 + t51 * t80 - t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 * t79, t103 * t58 + t71, 0, 0, 0, 0, 0, 0, -t57 + (t46 + t109) * qJD(4), -0.2e1 * t82 * qJD(4), -t43 - t139, t24 * t82 - t46 * t85 + t71, 0, 0, 0, 0, 0, 0, t89 - t130, -t126 + t144, (t14 - t131) * t77 + t34 * t102 + t117, t140 * t77 - t19 * t46 + t95 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, -t43 + t139, 0, -t122, t57 + (t46 - t109) * qJD(4), 0, -t55 * t46 - t88, (qJD(3) + t55) * t82, 0, 0, t101 * t34 - t133, (-t14 - t131) * t77 - t145 * t125 + t117, -t126 - t144, t102 * t32 - t132, t89 + t130, -t145 * t46, -pkin(4) * t15 - t12 * t77 - t24 * t32 + t93 * t46 + (-pkin(7) * t113 - t7) * t145 + t84 * t76, pkin(4) * t14 + t134 - t24 * t34 + t6 * t46 + (pkin(7) * t114 + t8) * t145 + t84 * t77, t32 * t8 + t34 * t7 + ((-t15 + t115) * pkin(7) + t95) * t77 + ((qJD(5) * t32 - t14) * pkin(7) - t140) * t76, -t12 * pkin(4) + pkin(7) * t80 - t19 * t24 - t6 * t8 + t7 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, -t32 ^ 2 + t34 ^ 2, t145 * t32 - t14, -t127, t145 * t34 - t15, t38, -t19 * t34 + t140, t19 * t32 - t95, 0, 0;];
tauc_reg = t5;
