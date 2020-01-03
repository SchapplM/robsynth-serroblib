% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:23
% EndTime: 2019-12-31 17:08:25
% DurationCPUTime: 0.68s
% Computational Cost: add. (794->138), mult. (1951->212), div. (0->0), fcn. (1094->4), ass. (0->100)
t72 = sin(qJ(4));
t100 = qJD(4) * t72;
t75 = cos(qJ(2));
t101 = qJD(2) * t75;
t73 = sin(qJ(2));
t74 = cos(qJ(4));
t99 = qJD(4) * t74;
t127 = -t75 * t100 + t72 * t101 + t73 * t99;
t37 = t73 * t72 + t75 * t74;
t29 = t37 * qJD(1);
t67 = qJD(2) - qJD(4);
t82 = t37 * qJD(4);
t97 = qJD(1) * qJD(2);
t93 = t75 * t97;
t94 = t73 * t97;
t7 = qJD(1) * t82 - t72 * t94 - t74 * t93;
t124 = t29 * t67 + t7;
t104 = qJD(1) * t73;
t64 = pkin(5) * t104;
t126 = -pkin(6) * t104 + qJD(3) + t64;
t102 = qJD(2) * t73;
t120 = pkin(5) - pkin(6);
t41 = t120 * t102;
t68 = qJD(2) * qJD(3);
t24 = -qJD(1) * t41 + t68;
t63 = pkin(5) * t93;
t34 = -pkin(6) * t93 + t63;
t121 = pkin(2) + pkin(3);
t21 = -t121 * qJD(2) + t126;
t103 = qJD(1) * t75;
t65 = pkin(5) * t103;
t42 = -pkin(6) * t103 + t65;
t69 = qJD(2) * qJ(3);
t32 = t42 + t69;
t6 = t72 * t21 + t74 * t32;
t2 = -qJD(4) * t6 - t72 * t24 + t74 * t34;
t125 = -t6 * t67 + t2;
t31 = -t72 * t103 + t74 * t104;
t8 = t127 * qJD(1) - t74 * t94;
t123 = t31 * t67 + t8;
t5 = t74 * t21 - t72 * t32;
t119 = t5 * t67;
t46 = t74 * qJ(3) - t121 * t72;
t117 = -t46 * qJD(4) - t126 * t72 - t74 * t42;
t92 = t73 * qJ(3) + pkin(1);
t35 = t121 * t75 + t92;
t20 = t35 * qJD(1);
t116 = t20 * t31;
t114 = t31 * t29;
t78 = qJD(1) ^ 2;
t112 = t75 * t78;
t77 = qJD(2) ^ 2;
t111 = t77 * t73;
t66 = t77 * t75;
t45 = -t72 * qJ(3) - t121 * t74;
t110 = t45 * qJD(4) + t126 * t74 - t72 * t42;
t70 = t73 ^ 2;
t107 = t75 ^ 2 - t70;
t106 = qJ(3) * t75;
t105 = qJD(2) * pkin(2);
t85 = pkin(2) * t73 - t106;
t98 = t73 * qJD(3);
t28 = t85 * qJD(2) - t98;
t19 = qJD(1) * t28;
t48 = -t75 * pkin(2) - t92;
t33 = qJD(1) * t48;
t54 = t120 * t75;
t95 = t29 ^ 2 - t31 ^ 2;
t90 = 0.2e1 * t33;
t89 = -0.2e1 * pkin(1) * t97;
t88 = qJD(3) - t105;
t87 = t67 ^ 2;
t86 = t73 * t93;
t53 = t120 * t73;
t14 = t74 * t53 - t72 * t54;
t15 = t72 * t53 + t74 * t54;
t84 = -pkin(5) * t77 - 0.2e1 * t19;
t83 = -t121 * t73 + t106;
t1 = -t32 * t100 + t21 * t99 + t74 * t24 + t72 * t34;
t80 = -t20 * t29 + t1;
t17 = t83 * qJD(2) + t98;
t44 = -pkin(5) * t94 + t68;
t47 = t64 + t88;
t50 = t65 + t69;
t79 = t44 * t75 + (t47 * t75 + (-t50 + t65) * t73) * qJD(2);
t60 = t73 * t112;
t52 = -0.2e1 * t86;
t51 = 0.2e1 * t86;
t49 = t107 * t78;
t43 = qJD(2) * t54;
t39 = t85 * qJD(1);
t38 = -t75 * t72 + t73 * t74;
t36 = t107 * t97;
t27 = t83 * qJD(1);
t13 = t17 * qJD(1);
t12 = t37 * qJD(2) - t82;
t11 = -t74 * t102 + t127;
t4 = -t15 * qJD(4) + t72 * t41 + t74 * t43;
t3 = t14 * qJD(4) - t74 * t41 + t72 * t43;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0.2e1 * t36, t66, t52, -t111, 0, -pkin(5) * t66 + t73 * t89, pkin(5) * t111 + t75 * t89, 0, 0, t51, t66, -0.2e1 * t36, 0, t111, t52, t90 * t102 + t84 * t75, t79, -t90 * t101 + t84 * t73, t79 * pkin(5) + t19 * t48 + t33 * t28, t31 * t12 - t7 * t38, -t31 * t11 - t12 * t29 + t7 * t37 - t38 * t8, -t12 * t67, t29 * t11 + t8 * t37, t11 * t67, 0, t20 * t11 + t13 * t37 + t17 * t29 + t35 * t8 - t4 * t67, t20 * t12 + t13 * t38 + t17 * t31 + t3 * t67 - t35 * t7, -t1 * t37 - t6 * t11 - t5 * t12 + t14 * t7 - t15 * t8 - t2 * t38 - t3 * t29 - t4 * t31, t1 * t15 + t13 * t35 + t2 * t14 + t20 * t17 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t49, 0, t60, 0, 0, t78 * pkin(1) * t73, pkin(1) * t112, 0, 0, -t60, 0, t49, 0, 0, t60, (-t33 * t73 + t39 * t75) * qJD(1), ((t50 - t69) * t73 + (-t47 + t88) * t75) * qJD(1), 0.2e1 * t68 + (t33 * t75 + t39 * t73) * qJD(1), t44 * qJ(3) + t50 * qJD(3) - t33 * t39 + (t50 * t73 + (-t47 - t105) * t75) * qJD(1) * pkin(5), -t114, t95, t124, t114, t123, 0, -t117 * t67 - t27 * t29 + t116 - t2, t110 * t67 - t27 * t31 + t80, t45 * t7 - t46 * t8 + (-t6 - t117) * t31 + (t5 - t110) * t29, t1 * t46 + t110 * t6 + t117 * t5 + t2 * t45 - t20 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, -t70 * t78 - t77, -t50 * qJD(2) + t33 * t104 + t63, 0, 0, 0, 0, 0, 0, -t29 * t104 - t72 * t87, -t31 * t104 - t74 * t87, -t123 * t72 + t124 * t74, -t20 * t104 + t125 * t74 + (t1 + t119) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t95, -t124, -t114, -t123, 0, -t116 + t125, -t80 - t119, 0, 0;];
tauc_reg = t9;
