% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:27
% EndTime: 2019-12-31 22:22:31
% DurationCPUTime: 1.08s
% Computational Cost: add. (2088->149), mult. (4808->262), div. (0->0), fcn. (4611->8), ass. (0->104)
t67 = cos(qJ(4));
t106 = qJD(4) * t67;
t64 = sin(qJ(3));
t114 = t64 * t67;
t63 = sin(qJ(4));
t68 = cos(qJ(3));
t123 = ((t63 * t68 + t114) * qJD(3) + t64 * t106) * pkin(2);
t65 = sin(qJ(2));
t69 = cos(qJ(2));
t43 = t64 * t69 + t68 * t65;
t120 = pkin(6) + pkin(7);
t49 = t120 * t65;
t50 = t120 * t69;
t83 = -t68 * t49 - t64 * t50;
t24 = -t43 * pkin(8) + t83;
t92 = qJD(2) * t120;
t44 = t65 * t92;
t45 = t69 * t92;
t84 = -t68 * t44 - t64 * t45;
t122 = t24 * qJD(3) + t84;
t66 = cos(qJ(5));
t61 = t66 ^ 2;
t62 = sin(qJ(5));
t109 = t62 ^ 2 - t61;
t89 = t109 * qJD(5);
t103 = t69 * qJD(2);
t104 = t65 * qJD(2);
t42 = t64 * t65 - t68 * t69;
t82 = t64 * t49 - t68 * t50;
t25 = -t42 * pkin(8) - t82;
t20 = t63 * t24 + t67 * t25;
t22 = t82 * qJD(3) + t64 * t44 - t68 * t45;
t34 = (-qJD(2) - qJD(3)) * t42;
t71 = t34 * pkin(8) - t22;
t5 = t20 * qJD(4) + t67 * t71 + ((-t64 * t103 - t68 * t104) * pkin(8) + t122) * t63;
t3 = t5 * t62;
t19 = -t67 * t24 + t63 * t25;
t59 = qJD(5) * t66;
t118 = t19 * t59 + t3;
t32 = t67 * t42 + t63 * t43;
t75 = t43 * qJD(2);
t70 = -t43 * qJD(3) - t75;
t11 = -t32 * qJD(4) + t67 * t34 + t63 * t70;
t33 = -t63 * t42 + t67 * t43;
t117 = t33 * t11;
t116 = t33 * t66;
t12 = t33 * qJD(4) + t63 * t34 - t67 * t70;
t115 = t62 * t12;
t113 = t66 * t11;
t112 = t66 * t12;
t107 = qJD(4) * t63;
t56 = t68 * pkin(2) + pkin(3);
t29 = t56 * t107 + t123;
t102 = t63 * t64 * pkin(2);
t38 = -t67 * t56 - pkin(4) + t102;
t111 = t29 * t62 + t38 * t59;
t55 = -t67 * pkin(3) - pkin(4);
t96 = pkin(3) * t107;
t110 = t55 * t59 + t62 * t96;
t108 = pkin(2) * qJD(3);
t105 = qJD(5) * t62;
t101 = -0.2e1 * pkin(1) * qJD(2);
t100 = pkin(4) * t105;
t99 = pkin(4) * t59;
t58 = pkin(2) * t104;
t98 = t64 * t108;
t97 = t68 * t108;
t95 = pkin(3) * t106;
t93 = t62 * t59;
t57 = -t69 * pkin(2) - pkin(1);
t91 = -0.4e1 * t62 * t116;
t35 = t38 * t105;
t90 = -t29 * t66 + t35;
t37 = t42 * pkin(3) + t57;
t18 = t32 * pkin(4) - t33 * pkin(9) + t37;
t88 = t66 * t18 - t62 * t20;
t87 = t62 * t18 + t66 * t20;
t39 = pkin(2) * t114 + t63 * t56 + pkin(9);
t86 = t32 * t39 - t33 * t38;
t54 = t63 * pkin(3) + pkin(9);
t85 = t32 * t54 - t33 * t55;
t46 = t55 * t105;
t80 = -t66 * t96 + t46;
t79 = t62 * t11 + t33 * t59;
t78 = t33 * t105 - t113;
t77 = t32 * t105 - t112;
t76 = t57 * t43;
t28 = -t56 * t106 - t67 * t97 + (qJD(3) + qJD(4)) * t102;
t74 = t11 * t38 - t12 * t39 + t28 * t32 + t29 * t33;
t73 = t11 * t55 - t12 * t54 + (-t32 * t67 + t33 * t63) * qJD(4) * pkin(3);
t27 = -t70 * pkin(3) + t58;
t53 = 0.2e1 * t93;
t41 = -0.2e1 * t89;
t31 = t33 ^ 2;
t21 = -t83 * qJD(3) - t84;
t16 = t19 * t105;
t10 = t32 * t59 + t115;
t8 = t62 * t113 - t33 * t89;
t7 = t12 * pkin(4) - t11 * pkin(9) + t27;
t6 = qJD(5) * t91 - t109 * t11;
t4 = t25 * t107 + t63 * t71 - t67 * (-pkin(8) * t75 + t122) - t24 * t106;
t2 = -t87 * qJD(5) + t62 * t4 + t66 * t7;
t1 = -t88 * qJD(5) + t66 * t4 - t62 * t7;
t9 = [0, 0, 0, 0.2e1 * t65 * t103, 0.2e1 * (-t65 ^ 2 + t69 ^ 2) * qJD(2), 0, 0, 0, t65 * t101, t69 * t101, 0.2e1 * t43 * t34, -0.2e1 * t34 * t42 + 0.2e1 * t43 * t70, 0, 0, 0, 0.2e1 * qJD(3) * t76 + 0.2e1 * (t65 * pkin(2) * t42 + t76) * qJD(2), 0.2e1 * t57 * t34 + 0.2e1 * t43 * t58, 0.2e1 * t117, -0.2e1 * t11 * t32 - 0.2e1 * t33 * t12, 0, 0, 0, 0.2e1 * t37 * t12 + 0.2e1 * t27 * t32, 0.2e1 * t37 * t11 + 0.2e1 * t27 * t33, 0.2e1 * t61 * t117 - 0.2e1 * t31 * t93, t11 * t91 + 0.2e1 * t31 * t89, 0.2e1 * t33 * t112 - 0.2e1 * t78 * t32, -0.2e1 * t33 * t115 - 0.2e1 * t79 * t32, 0.2e1 * t32 * t12, 0.2e1 * t88 * t12 + 0.2e1 * t79 * t19 + 0.2e1 * t2 * t32 + 0.2e1 * t33 * t3, 0.2e1 * t1 * t32 + 0.2e1 * t5 * t116 - 0.2e1 * t87 * t12 - 0.2e1 * t78 * t19; 0, 0, 0, 0, 0, t103, -t104, 0, -pkin(6) * t103, pkin(6) * t104, 0, 0, t34, t70, 0, t22, t21, 0, 0, t11, -t12, 0, -t5, t4, t8, t6, t10, -t77, 0, t16 + (-t86 * qJD(5) - t5) * t66 + t74 * t62, t86 * t105 + t74 * t66 + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t98, -0.2e1 * t97, 0, 0, 0, 0, 0, -0.2e1 * t29, 0.2e1 * t28, t53, t41, 0, 0, 0, 0.2e1 * t90, 0.2e1 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t70, 0, t22, t21, 0, 0, t11, -t12, 0, -t5, t4, t8, t6, t10, -t77, 0, t16 + (-t85 * qJD(5) - t5) * t66 + t73 * t62, t85 * t105 + t73 * t66 + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t97, 0, 0, 0, 0, 0, (-pkin(3) - t56) * t107 - t123, t28 - t95, t53, t41, 0, 0, 0, t35 + t46 + (-t29 - t96) * t66, t110 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t96, -0.2e1 * t95, t53, t41, 0, 0, 0, 0.2e1 * t80, 0.2e1 * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t12, 0, -t5, t4, t8, t6, t10, -t77, 0, t16 + (-pkin(4) * t11 - pkin(9) * t12) * t62 + (-t5 + (-pkin(4) * t33 - pkin(9) * t32) * qJD(5)) * t66, t78 * pkin(4) + t77 * pkin(9) + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t28, t53, t41, 0, 0, 0, t90 - t100, -t99 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t95, t53, t41, 0, 0, 0, t80 - t100, -t99 + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t41, 0, 0, 0, -0.2e1 * t100, -0.2e1 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t79, t12, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t105, 0, t62 * t28 - t39 * t59, t39 * t105 + t66 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t105, 0, -t54 * t59 - t62 * t95, t54 * t105 - t66 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t105, 0, -pkin(9) * t59, pkin(9) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
