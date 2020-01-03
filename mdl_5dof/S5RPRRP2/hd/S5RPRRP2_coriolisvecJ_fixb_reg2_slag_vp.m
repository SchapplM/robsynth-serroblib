% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP2
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:29
% EndTime: 2020-01-03 11:45:32
% DurationCPUTime: 0.69s
% Computational Cost: add. (1131->156), mult. (2398->202), div. (0->0), fcn. (1288->6), ass. (0->101)
t59 = cos(pkin(8)) * pkin(1) + pkin(2);
t50 = t59 * qJD(1);
t105 = qJD(3) * t50;
t73 = sin(qJ(3));
t75 = cos(qJ(3));
t125 = sin(pkin(8)) * pkin(1);
t98 = qJD(1) * t125;
t86 = qJD(3) * t98;
t32 = t75 * t105 - t73 * t86;
t127 = -qJD(2) * qJD(4) - t32;
t72 = sin(qJ(4));
t102 = t72 * qJD(2);
t74 = cos(qJ(4));
t36 = t73 * t50 + t75 * t98;
t67 = qJD(1) + qJD(3);
t26 = t67 * pkin(7) + t36;
t90 = qJ(5) * t67 + t26;
t81 = t90 * t74;
t11 = t81 + t102;
t63 = t74 * qJD(2);
t16 = -t72 * t26 + t63;
t17 = t74 * t26 + t102;
t6 = -t17 * qJD(4) - t72 * t32;
t104 = qJD(4) * t72;
t99 = t26 * t104 + t127 * t74;
t126 = -t99 * t74 + (-t16 * t74 - t17 * t72) * qJD(4) - t6 * t72;
t110 = t75 * t125 + t73 * t59;
t91 = -t73 * t125 + t75 * t59;
t68 = t72 ^ 2;
t124 = pkin(4) * t68;
t123 = t74 * pkin(4);
t10 = -t72 * t90 + t63;
t107 = qJD(4) * pkin(4);
t9 = t10 + t107;
t122 = -t10 + t9;
t37 = t91 * qJD(3);
t121 = t37 * t67;
t38 = t110 * qJD(3);
t120 = t38 * t67;
t66 = t67 ^ 2;
t119 = t66 * t74;
t118 = t67 * t72;
t117 = t67 * t74;
t76 = qJD(4) ^ 2;
t115 = t76 * t72;
t114 = -qJ(5) - pkin(7);
t103 = qJD(4) * t74;
t35 = t75 * t50 - t73 * t98;
t61 = -pkin(3) - t123;
t15 = t61 * t67 + qJD(5) - t35;
t33 = t73 * t105 + t75 * t86;
t96 = t67 * t104;
t18 = pkin(4) * t96 + t33;
t113 = t15 * t103 + t18 * t72;
t25 = -t67 * pkin(3) - t35;
t112 = t25 * t103 + t33 * t72;
t111 = t35 * t104 + t36 * t117;
t69 = t74 ^ 2;
t109 = t68 - t69;
t108 = t68 + t69;
t41 = pkin(7) + t110;
t106 = -qJ(5) - t41;
t62 = t74 * qJD(5);
t101 = -qJD(5) - t15;
t97 = pkin(4) * t104;
t95 = t67 * t103;
t93 = qJ(5) * t104;
t2 = (-t93 + t62) * t67 - t99;
t3 = (-qJD(5) * t67 - t32) * t72 - t11 * qJD(4);
t94 = t2 * t74 - t3 * t72;
t92 = t108 * t35;
t89 = qJD(4) * t114;
t88 = qJD(4) * t106;
t87 = t72 * t95;
t40 = -pkin(3) - t91;
t85 = t11 * t74 - t72 * t9;
t84 = -t11 * t72 - t74 * t9;
t83 = t16 * t72 - t17 * t74;
t82 = t41 * t76 + t120;
t80 = qJD(4) * (t40 * t67 - t37);
t65 = t74 * qJ(5);
t64 = t76 * t74;
t53 = t72 * t119;
t52 = t74 * pkin(7) + t65;
t51 = t114 * t72;
t47 = -0.2e1 * t87;
t46 = 0.2e1 * t87;
t44 = t109 * t66;
t43 = -t72 * qJD(5) + t74 * t89;
t42 = t72 * t89 + t62;
t39 = -0.2e1 * t109 * t67 * qJD(4);
t34 = t40 - t123;
t31 = t38 + t97;
t30 = t74 * t41 + t65;
t29 = t106 * t72;
t28 = t35 * t103;
t19 = t25 * t104;
t12 = t15 * t104;
t8 = (-qJD(5) - t37) * t72 + t74 * t88;
t7 = t74 * t37 + t72 * t88 + t62;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33 - t120, -t32 - t121, 0, t32 * t110 - t33 * t91 - t35 * t38 + t36 * t37, t46, t39, t64, t47, -t115, 0, t19 + t72 * t80 + (-t33 - t82) * t74, t72 * t82 + t74 * t80 + t112, t108 * t121 + t126, t126 * t41 + t25 * t38 + t33 * t40 - t37 * t83, t46, t39, t64, t47, -t115, 0, t12 + (-t31 * t67 - t18) * t74 + (t34 * t118 + t8) * qJD(4), t31 * t118 + (t34 * t117 - t7) * qJD(4) + t113, (t7 * t74 - t72 * t8) * t67 + ((-t29 * t74 - t30 * t72) * t67 + t84) * qJD(4) + t94, t11 * t7 + t15 * t31 + t18 * t34 + t2 * t30 + t3 * t29 + t9 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t64, 0, -qJD(4) * t83 + t6 * t74 - t72 * t99, 0, 0, 0, 0, 0, 0, -t115, -t64, 0, qJD(4) * t85 + t2 * t72 + t3 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t67 - t33, t35 * t67 - t32, 0, 0, t46, t39, t64, t47, -t115, 0, -pkin(3) * t96 + t19 + (-pkin(7) * t76 - t33) * t74 + t111, pkin(7) * t115 + t28 + (-pkin(3) * t103 - t36 * t72) * t67 + t112, -t67 * t92 + t126, -t33 * pkin(3) + pkin(7) * t126 - t25 * t36 + t35 * t83, t46, t39, t64, t47, -t115, 0, -t18 * t74 + t12 + (t43 + (t61 - t123) * t118) * qJD(4) + t111, -t36 * t118 + t28 + (-t42 + (t61 * t74 + t124) * t67) * qJD(4) + t113, t84 * qJD(4) + (t42 * t74 - t43 * t72 - t92 + (-t51 * t74 - t52 * t72) * qJD(4)) * t67 + t94, t18 * t61 + t2 * t52 + t3 * t51 + (t35 * t72 + t43) * t9 + (-t36 + t97) * t15 + (-t35 * t74 + t42) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t44, 0, t53, 0, 0, (-t25 * t67 - t32) * t72, t16 * qJD(4) - t25 * t117 + t99, 0, 0, -t53, t44, 0, t53, 0, 0, (t11 - t81) * qJD(4) + (pkin(4) * t119 + t101 * t67 + t127) * t72, -t66 * t124 + t10 * qJD(4) + (t101 * t74 + t93) * t67 + t99, (-t107 + t122) * t117, t122 * t11 + (-t15 * t118 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t96, 0.2e1 * t95, -t108 * t66, -t67 * t85 + t18;];
tauc_reg = t1;
