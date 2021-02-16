% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP12
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
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:25:53
% EndTime: 2021-01-15 19:25:58
% DurationCPUTime: 0.79s
% Computational Cost: add. (552->137), mult. (1262->239), div. (0->0), fcn. (880->4), ass. (0->91)
t40 = cos(qJ(3));
t38 = sin(qJ(3));
t77 = t38 * qJD(3);
t102 = qJ(5) * t77 - t40 * qJD(5);
t101 = 2 * qJD(2);
t100 = 2 * qJD(4);
t39 = cos(qJ(4));
t99 = t39 * pkin(4);
t98 = t40 * pkin(7);
t37 = sin(qJ(4));
t79 = qJD(4) * t40;
t67 = t39 * t79;
t17 = -t37 * t77 + t67;
t41 = -pkin(1) - pkin(6);
t64 = t41 * t77;
t7 = t17 * pkin(4) + t64;
t97 = t7 * t37;
t96 = t7 * t39;
t89 = -qJ(5) - pkin(7);
t23 = t89 * t37;
t95 = t23 * t40;
t24 = t89 * t39;
t94 = t24 * t40;
t30 = -pkin(3) - t99;
t93 = t30 * t39;
t92 = t37 * t41;
t91 = t38 * t41;
t90 = t40 * t41;
t51 = t38 * pkin(3) - t98;
t47 = qJ(2) + t51;
t19 = t37 * t47;
t28 = t39 * t91;
t88 = t28 + t19;
t33 = t37 ^ 2;
t35 = t39 ^ 2;
t87 = t33 - t35;
t86 = t33 + t35;
t34 = t38 ^ 2;
t36 = t40 ^ 2;
t85 = t34 - t36;
t84 = t34 + t36;
t83 = qJ(5) * t40;
t21 = (pkin(4) * t37 - t41) * t40;
t82 = qJD(3) * t21;
t81 = qJD(3) * t39;
t80 = qJD(4) * t37;
t32 = qJD(4) * t39;
t78 = qJD(4) * t41;
t76 = t40 * qJD(3);
t74 = qJ(2) * qJD(3);
t73 = -0.2e1 * pkin(3) * qJD(4);
t20 = t39 * t47;
t52 = pkin(3) * t40 + pkin(7) * t38;
t45 = t52 * qJD(3) + qJD(2);
t61 = t41 * t76;
t72 = -qJD(4) * t20 - t37 * t45 - t39 * t61;
t71 = pkin(4) * t80;
t70 = t39 * t83;
t69 = t37 * t79;
t68 = t37 * t78;
t66 = t21 * t80;
t65 = t37 * t32;
t63 = t39 * t77;
t62 = t38 * t76;
t59 = -t30 + t99;
t58 = pkin(4) - t92;
t57 = t87 * qJD(4);
t56 = t85 * qJD(3);
t55 = 0.2e1 * t62;
t54 = t37 * t61;
t53 = t37 * t63;
t5 = t58 * t38 + t20 - t70;
t6 = -t37 * t83 + t88;
t50 = t37 * t6 + t39 * t5;
t49 = t37 * t5 - t39 * t6;
t48 = pkin(4) * t33 + t93;
t46 = t63 + t69;
t16 = t38 * t32 + t37 * t76;
t8 = -t39 * qJD(5) - t89 * t80;
t9 = -t37 * qJD(5) + t89 * t32;
t44 = -t9 * t37 - t8 * t39 + (-t23 * t39 + t24 * t37) * qJD(4);
t43 = -t88 * qJD(4) + t39 * t45;
t42 = qJ(5) * t69 + t102 * t39 + t43;
t15 = t84 * t32;
t13 = t38 * t80 - t39 * t76;
t12 = t84 * t80;
t4 = t43 - t54;
t3 = t38 * t68 + t72;
t2 = t70 * qJD(4) + t72 + (t91 * qJD(4) - t102) * t37;
t1 = t58 * t76 + t42;
t10 = [0, 0, 0, 0, t101, qJ(2) * t101, -0.2e1 * t62, 0.2e1 * t56, 0, 0, 0, 0.2e1 * qJD(2) * t38 + 0.2e1 * t40 * t74, 0.2e1 * qJD(2) * t40 - 0.2e1 * t38 * t74, -0.2e1 * t35 * t62 - 0.2e1 * t36 * t65, t87 * t36 * t100 + 0.4e1 * t40 * t53, -0.2e1 * t38 * t69 - 0.2e1 * t85 * t81, 0.2e1 * t37 * t56 - 0.2e1 * t38 * t67, t55, -0.2e1 * t36 * t39 * t78 + 0.2e1 * t20 * t76 + 0.2e1 * (t4 + t54) * t38, 0.2e1 * t36 * t68 + 0.2e1 * t3 * t38 + 0.2e1 * (-t19 + t28) * t76, 0.2e1 * (-t37 * t82 + t1) * t38 + 0.2e1 * (qJD(3) * t5 + t21 * t32 + t97) * t40, 0.2e1 * (-t21 * t81 + t2) * t38 + 0.2e1 * (-qJD(3) * t6 - t66 + t96) * t40, 0.2e1 * t50 * t77 + 0.2e1 * (t49 * qJD(4) - t1 * t39 + t2 * t37) * t40, 0.2e1 * t5 * t1 - 0.2e1 * t6 * t2 + 0.2e1 * t21 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t12, -t15, t12, 0, (-t49 * qJD(3) - t7) * t40 + (-t50 * qJD(4) - t1 * t37 - t2 * t39 + t82) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t86) * t55; 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76, 0, -t64, -t61, -t40 * t57 - t53, -0.4e1 * t40 * t65 + t87 * t77, t16, -t13, 0, (-t37 * t90 - t52 * t39) * qJD(4) + (t51 * t37 - t28) * qJD(3), (t52 * t37 - t39 * t90) * qJD(4) + (-t39 * t98 + (pkin(3) * t39 + t92) * t38) * qJD(3), t9 * t38 - t96 + (-t30 * t37 * t38 + t95) * qJD(3) + (t21 * t37 + t48 * t40) * qJD(4), t97 + t8 * t38 + (-t38 * t93 + t94) * qJD(3) + (t59 * t40 * t37 + t21 * t39) * qJD(4), (t23 * t77 - t40 * t9 - t2 + (-t5 + t94) * qJD(4)) * t39 + (-t24 * t77 + t40 * t8 - t1 + (-t6 + t95) * qJD(4)) * t37, pkin(4) * t66 + t1 * t23 + t2 * t24 + t7 * t30 + t5 * t9 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76, 0, 0, 0, 0, 0, -t46, -t17, -t46, -t17, t86 * t76, (-t71 + (-t23 * t37 - t24 * t39) * qJD(3)) * t40 + (qJD(3) * t30 + t44) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t65, -0.2e1 * t57, 0, 0, 0, t37 * t73, t39 * t73, -0.2e1 * t59 * t80, t48 * t100, 0.2e1 * t44, 0.2e1 * t23 * t9 + 0.2e1 * t24 * t8 + 0.2e1 * t30 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t17, t76, t4, t3, (0.2e1 * pkin(4) - t92) * t76 + t42, t2, t46 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t13, -t16, t13, 0, -t16 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t80, 0, -pkin(7) * t32, pkin(7) * t80, t9, t8, -pkin(4) * t32, t9 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t46, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t32, 0, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
