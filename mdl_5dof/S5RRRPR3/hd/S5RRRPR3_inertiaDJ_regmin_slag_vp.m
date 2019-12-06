% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:58
% EndTime: 2019-12-05 18:43:00
% DurationCPUTime: 0.52s
% Computational Cost: add. (868->101), mult. (2011->169), div. (0->0), fcn. (1787->8), ass. (0->92)
t94 = sin(pkin(9));
t125 = pkin(3) * t94;
t97 = sin(qJ(3));
t112 = t97 * qJD(3);
t100 = cos(qJ(3));
t92 = t100 * qJD(3);
t95 = cos(pkin(9));
t66 = -t94 * t112 + t95 * t92;
t124 = t66 * pkin(8);
t70 = t95 * t100 - t94 * t97;
t123 = t70 * pkin(4);
t71 = t94 * t100 + t95 * t97;
t122 = t71 * pkin(8);
t96 = sin(qJ(5));
t99 = cos(qJ(5));
t103 = t99 * t70 - t96 * t71;
t41 = t96 * t70 + t99 * t71;
t65 = t71 * qJD(3);
t16 = t41 * qJD(5) + t99 * t65 + t96 * t66;
t89 = pkin(3) * t112;
t51 = t65 * pkin(4) + t89;
t114 = pkin(1) * qJD(2);
t98 = sin(qJ(2));
t90 = t98 * t114;
t47 = t51 + t90;
t101 = cos(qJ(2));
t118 = t101 * pkin(1);
t88 = -t100 * pkin(3) - pkin(2);
t79 = t88 - t118;
t50 = t79 - t123;
t121 = -t103 * t47 + t50 * t16;
t15 = t103 * qJD(5) - t96 * t65 + t99 * t66;
t120 = t50 * t15 + t47 * t41;
t52 = t88 - t123;
t119 = t52 * t15 + t51 * t41;
t117 = -qJ(4) - pkin(7);
t116 = -t103 * t51 + t52 * t16;
t109 = t101 * t114;
t104 = t100 * t109;
t86 = t98 * pkin(1) + pkin(7);
t113 = -qJ(4) - t86;
t105 = qJD(3) * t113;
t91 = t100 * qJD(4);
t45 = t97 * t105 + t104 + t91;
t46 = (-qJD(4) - t109) * t97 + t100 * t105;
t18 = t95 * t45 + t94 * t46;
t108 = qJD(3) * t117;
t63 = t97 * t108 + t91;
t64 = -t97 * qJD(4) + t100 * t108;
t35 = t95 * t63 + t94 * t64;
t68 = t113 * t97;
t93 = t100 * qJ(4);
t69 = t100 * t86 + t93;
t39 = t94 * t68 + t95 * t69;
t80 = t117 * t97;
t81 = t100 * pkin(7) + t93;
t49 = t94 * t80 + t95 * t81;
t87 = -pkin(2) - t118;
t115 = t87 * t92 + t97 * t90;
t111 = pkin(2) * t112;
t110 = pkin(2) * t92;
t17 = -t94 * t45 + t95 * t46;
t34 = -t94 * t63 + t95 * t64;
t38 = t95 * t68 - t94 * t69;
t48 = t95 * t80 - t94 * t81;
t107 = -t17 * t71 + t18 * t70 - t38 * t66 - t39 * t65;
t106 = -t34 * t71 + t35 * t70 - t48 * t66 - t49 * t65;
t102 = -t100 * t90 + t87 * t112;
t85 = t95 * pkin(3) + pkin(4);
t83 = 0.2e1 * t97 * t92;
t76 = t90 + t89;
t75 = 0.2e1 * (t100 ^ 2 - t97 ^ 2) * qJD(3);
t67 = t70 * pkin(8);
t62 = t65 * pkin(8);
t57 = (-t99 * t125 - t85 * t96) * qJD(5);
t56 = (t96 * t125 - t85 * t99) * qJD(5);
t37 = t67 + t49;
t36 = t48 - t122;
t33 = (-t65 * t94 - t66 * t95) * pkin(3);
t30 = t67 + t39;
t29 = t38 - t122;
t22 = -t62 + t35;
t21 = t34 - t124;
t12 = -t62 + t18;
t11 = t17 - t124;
t6 = 0.2e1 * t41 * t15;
t5 = t99 * t21 - t96 * t22 + (-t36 * t96 - t37 * t99) * qJD(5);
t4 = -t96 * t21 - t99 * t22 + (-t36 * t99 + t37 * t96) * qJD(5);
t3 = 0.2e1 * t103 * t15 - 0.2e1 * t41 * t16;
t2 = t99 * t11 - t96 * t12 + (-t29 * t96 - t30 * t99) * qJD(5);
t1 = -t96 * t11 - t99 * t12 + (-t29 * t99 + t30 * t96) * qJD(5);
t7 = [0, 0, 0, 0, -0.2e1 * t90, -0.2e1 * t109, t83, t75, 0, 0, 0, 0.2e1 * t102, 0.2e1 * t115, 0.2e1 * t107, 0.2e1 * t38 * t17 + 0.2e1 * t39 * t18 + 0.2e1 * t79 * t76, t6, t3, 0, 0, 0, 0.2e1 * t121, 0.2e1 * t120; 0, 0, 0, 0, -t90, -t109, t83, t75, 0, 0, 0, t102 - t111, -t110 + t115, t106 + t107, t17 * t48 + t18 * t49 + t38 * t34 + t39 * t35 + t76 * t88 + t79 * t89, t6, t3, 0, 0, 0, t116 + t121, t119 + t120; 0, 0, 0, 0, 0, 0, t83, t75, 0, 0, 0, -0.2e1 * t111, -0.2e1 * t110, 0.2e1 * t106, 0.2e1 * t48 * t34 + 0.2e1 * t49 * t35 + 0.2e1 * t88 * t89, t6, t3, 0, 0, 0, 0.2e1 * t116, 0.2e1 * t119; 0, 0, 0, 0, 0, 0, 0, 0, t92, -t112, 0, -t97 * t109 - t86 * t92, t86 * t112 - t104, t33, (t17 * t95 + t18 * t94) * pkin(3), 0, 0, t15, -t16, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, t92, -t112, 0, -pkin(7) * t92, pkin(7) * t112, t33, (t34 * t95 + t35 * t94) * pkin(3), 0, 0, t15, -t16, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t57, 0.2e1 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, 0, 0, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
