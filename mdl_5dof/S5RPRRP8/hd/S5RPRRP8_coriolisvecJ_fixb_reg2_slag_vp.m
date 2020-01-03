% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:30
% EndTime: 2019-12-31 18:47:32
% DurationCPUTime: 0.68s
% Computational Cost: add. (1206->133), mult. (1971->164), div. (0->0), fcn. (845->4), ass. (0->97)
t56 = -pkin(1) - pkin(2);
t40 = t56 * qJD(1) + qJD(2);
t53 = sin(qJ(3));
t103 = t53 * t40;
t55 = cos(qJ(3));
t88 = qJD(1) * qJ(2);
t28 = t55 * t88 + t103;
t85 = qJD(1) - qJD(3);
t110 = t28 * t85;
t89 = t53 * qJD(2);
t94 = qJD(3) * t55;
t15 = (qJ(2) * t94 + t89) * qJD(1) + qJD(3) * t103;
t121 = -t110 - t15;
t120 = t85 ^ 2;
t52 = sin(qJ(4));
t50 = t52 ^ 2;
t54 = cos(qJ(4));
t51 = t54 ^ 2;
t95 = t50 + t51;
t119 = qJD(4) * t85;
t42 = t53 * t88;
t86 = qJD(1) * qJD(2);
t14 = -qJD(3) * t42 + t40 * t94 + t55 * t86;
t27 = t55 * t40 - t42;
t118 = t85 * t27 + t14;
t117 = 0.2e1 * t55 * t119;
t77 = t95 * t14;
t97 = t55 * qJ(2) + t53 * t56;
t26 = t97 * qJD(3) + t89;
t116 = t26 * t85 + t15;
t75 = -t53 * qJ(2) + t55 * t56;
t66 = pkin(4) * t52 - qJ(5) * t54;
t30 = t66 * qJD(4) - t52 * qJD(5);
t20 = -pkin(7) * t85 + t28;
t104 = t52 * t20;
t12 = t54 * t14;
t1 = t12 + (qJD(5) - t104) * qJD(4);
t102 = t54 * t20;
t87 = qJD(4) * qJ(5);
t10 = t87 + t102;
t11 = t52 * t14;
t92 = qJD(4) * t54;
t2 = t20 * t92 + t11;
t74 = -qJD(4) * pkin(4) + qJD(5);
t9 = t74 + t104;
t59 = -(t10 * t52 - t54 * t9) * qJD(4) + t1 * t54 + t2 * t52;
t96 = -t50 + t51;
t115 = -0.2e1 * t96 * t119;
t57 = qJD(4) ^ 2;
t114 = pkin(7) * t57;
t113 = t85 * pkin(3);
t25 = t55 * qJD(2) + t75 * qJD(3);
t112 = t25 * t85;
t109 = t30 * t85;
t35 = -pkin(7) + t97;
t108 = t35 * t57;
t38 = -t54 * pkin(4) - t52 * qJ(5) - pkin(3);
t107 = t38 * t85;
t106 = t85 * t52;
t105 = t85 * t54;
t100 = t95 * t112;
t93 = qJD(4) * t52;
t99 = -t28 * t105 + t27 * t93;
t98 = t28 - t30;
t83 = 0.2e1 * t86;
t5 = t15 - t109;
t82 = -t5 - t114;
t19 = -t27 + t113;
t81 = t19 + t113;
t80 = -t15 - t114;
t7 = -t27 - t107;
t79 = t7 - t107;
t78 = t95 * t27;
t72 = t92 * t106;
t22 = -t38 - t75;
t71 = -t22 * t85 - t25 - t7;
t70 = t85 * t78;
t68 = t10 * t54 + t52 * t9;
t8 = t26 - t30;
t65 = t8 * t85 - t108 + t5;
t64 = t108 - t116;
t34 = pkin(3) - t75;
t63 = qJD(4) * (-t34 * t85 - t19 - t25);
t61 = (-t57 - t120) * t53;
t60 = t85 * t95;
t58 = qJD(1) ^ 2;
t45 = t57 * t54;
t44 = t57 * t52;
t39 = t52 * t120 * t54;
t37 = 0.2e1 * t72;
t36 = -0.2e1 * t72;
t32 = t96 * t120;
t31 = t66 * t85;
t6 = t60 * t55 * t85;
t4 = t52 * t117 + t54 * t61;
t3 = -t54 * t117 + t52 * t61;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, qJ(2) * t83, 0, 0, 0, 0, 0, 0, t116, t14 + t112, 0, t14 * t97 - t15 * t75 + t28 * t25 - t27 * t26, t37, -t115, -t45, t36, t44, 0, t52 * t63 - t54 * t64, t52 * t64 + t54 * t63, -t100 - t77, t95 * t25 * t20 + t15 * t34 + t19 * t26 + t35 * t77, t37, -t45, t115, 0, -t44, t36, t65 * t54 + t71 * t93, -t100 - t59, t65 * t52 - t71 * t92, t5 * t22 + t25 * t68 + t35 * t59 + t7 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t58 * qJ(2), 0, 0, 0, 0, 0, 0, -t53 * t120, -t55 * t120, 0, t118 * t53 + t121 * t55, 0, 0, 0, 0, 0, 0, t4, -t3, t6, (-t85 * t19 + t77) * t53 + (-t20 * t60 - t15) * t55, 0, 0, 0, 0, 0, 0, t4, t6, t3, (-t85 * t68 - t5) * t55 + (-t85 * t7 + t59) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, -t118, 0, 0, t36, t115, t45, t37, -t44, 0, t80 * t54 + t81 * t93 + t99, (-t80 + t110) * t52 + (t27 + t81) * t92, t70 + t77, -t15 * pkin(3) + pkin(7) * t77 - t19 * t28 - t20 * t78, t36, t45, -t115, 0, t44, t37, t79 * t93 + (t82 + t109) * t54 + t99, t70 + t59, (-t27 - t79) * t92 + (-t85 * t98 + t82) * t52, t59 * pkin(7) - t68 * t27 + t5 * t38 - t98 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t32, 0, t39, 0, 0, t19 * t106 - t11, t19 * t105 - t12, 0, 0, -t39, 0, t32, 0, 0, t39, -t11 - (-t31 * t54 - t52 * t7) * t85, -((t10 - t87) * t52 + (t74 - t9) * t54) * t85, 0.2e1 * qJD(4) * qJD(5) + t12 - (-t31 * t52 + t54 * t7) * t85, -t9 * t102 - t2 * pkin(4) + t1 * qJ(5) + t7 * t31 + (qJD(5) + t104) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, -t120 * t50 - t57, -t10 * qJD(4) - t7 * t106 + t2;];
tauc_reg = t13;
