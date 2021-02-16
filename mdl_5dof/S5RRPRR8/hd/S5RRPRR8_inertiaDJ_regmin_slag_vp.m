% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:11
% EndTime: 2021-01-15 21:35:16
% DurationCPUTime: 0.67s
% Computational Cost: add. (1385->116), mult. (3184->224), div. (0->0), fcn. (3144->8), ass. (0->86)
t105 = cos(qJ(4));
t63 = sin(pkin(9));
t107 = pkin(2) * t63;
t66 = sin(qJ(4));
t64 = cos(pkin(9));
t87 = t64 * pkin(2) + pkin(3);
t108 = -t105 * t87 + t66 * t107;
t67 = sin(qJ(2));
t98 = -qJ(3) - pkin(6);
t53 = t98 * t67;
t69 = cos(qJ(2));
t54 = t98 * t69;
t29 = t64 * t53 + t63 * t54;
t47 = t63 * t69 + t64 * t67;
t24 = -t47 * pkin(7) + t29;
t30 = t63 * t53 - t64 * t54;
t46 = t63 * t67 - t64 * t69;
t25 = -t46 * pkin(7) + t30;
t14 = -t105 * t24 + t66 * t25;
t68 = cos(qJ(5));
t62 = t68 ^ 2;
t65 = sin(qJ(5));
t96 = t65 ^ 2 - t62;
t82 = t96 * qJD(5);
t15 = t105 * t25 + t66 * t24;
t83 = qJD(2) * t98;
t41 = t69 * qJD(3) + t67 * t83;
t42 = -t67 * qJD(3) + t69 * t83;
t22 = -t63 * t41 + t64 * t42;
t93 = t69 * qJD(2);
t94 = t67 * qJD(2);
t45 = -t63 * t94 + t64 * t93;
t72 = t45 * pkin(7) - t22;
t23 = t64 * t41 + t63 * t42;
t44 = t47 * qJD(2);
t73 = -t44 * pkin(7) + t23;
t5 = t15 * qJD(4) + t105 * t72 + t66 * t73;
t3 = t5 * t65;
t60 = qJD(5) * t68;
t106 = t14 * t60 + t3;
t74 = -t105 * t46 - t66 * t47;
t16 = t74 * qJD(4) + t105 * t45 - t66 * t44;
t28 = t105 * t47 - t66 * t46;
t104 = t28 * t16;
t103 = t28 * t68;
t17 = t28 * qJD(4) + t105 * t44 + t66 * t45;
t102 = t65 * t17;
t100 = t68 * t16;
t99 = t68 * t17;
t71 = t105 * t107 + t66 * t87;
t37 = t71 * qJD(4);
t39 = -pkin(4) + t108;
t97 = t37 * t65 + t39 * t60;
t95 = qJD(5) * t65;
t91 = -0.2e1 * pkin(1) * qJD(2);
t90 = pkin(4) * t95;
t89 = pkin(4) * t60;
t59 = pkin(2) * t94;
t88 = t65 * t60;
t58 = -t69 * pkin(2) - pkin(1);
t31 = t44 * pkin(3) + t59;
t85 = -0.4e1 * t65 * t103;
t84 = -t37 * t68 + t39 * t95;
t35 = t46 * pkin(3) + t58;
t13 = -pkin(4) * t74 - t28 * pkin(8) + t35;
t81 = t68 * t13 - t65 * t15;
t80 = t65 * t13 + t68 * t15;
t40 = pkin(8) + t71;
t79 = -t28 * t39 - t40 * t74;
t77 = t65 * t16 + t28 * t60;
t76 = t28 * t95 - t100;
t10 = -t60 * t74 + t102;
t75 = -t74 * t95 - t99;
t36 = t108 * qJD(4);
t70 = t16 * t39 - t17 * t40 + t28 * t37 - t36 * t74;
t55 = 0.2e1 * t88;
t51 = -0.2e1 * t82;
t26 = t28 ^ 2;
t11 = t14 * t95;
t8 = t65 * t100 - t28 * t82;
t7 = t17 * pkin(4) - t16 * pkin(8) + t31;
t6 = qJD(5) * t85 - t96 * t16;
t4 = t14 * qJD(4) - t105 * t73 + t66 * t72;
t2 = -t80 * qJD(5) + t65 * t4 + t68 * t7;
t1 = -t81 * qJD(5) + t68 * t4 - t65 * t7;
t9 = [0, 0, 0, 0.2e1 * t67 * t93, 0.2e1 * (-t67 ^ 2 + t69 ^ 2) * qJD(2), 0, 0, 0, t67 * t91, t69 * t91, 0.2e1 * t58 * t44 + 0.2e1 * t46 * t59, 0.2e1 * t58 * t45 + 0.2e1 * t47 * t59, -0.2e1 * t22 * t47 - 0.2e1 * t23 * t46 - 0.2e1 * t29 * t45 - 0.2e1 * t30 * t44, 0.2e1 * t29 * t22 + 0.2e1 * t30 * t23 + 0.2e1 * t58 * t59, 0.2e1 * t104, 0.2e1 * t16 * t74 - 0.2e1 * t28 * t17, 0, 0, 0, 0.2e1 * t35 * t17 - 0.2e1 * t31 * t74, 0.2e1 * t35 * t16 + 0.2e1 * t31 * t28, 0.2e1 * t62 * t104 - 0.2e1 * t26 * t88, t16 * t85 + 0.2e1 * t26 * t82, 0.2e1 * t28 * t99 + 0.2e1 * t74 * t76, -0.2e1 * t28 * t102 + 0.2e1 * t74 * t77, -0.2e1 * t74 * t17, 0.2e1 * t77 * t14 + 0.2e1 * t81 * t17 - 0.2e1 * t2 * t74 + 0.2e1 * t28 * t3, -0.2e1 * t1 * t74 + 0.2e1 * t5 * t103 - 0.2e1 * t76 * t14 - 0.2e1 * t80 * t17; 0, 0, 0, 0, 0, t93, -t94, 0, -pkin(6) * t93, pkin(6) * t94, t22, -t23, (-t44 * t63 - t45 * t64) * pkin(2), (t22 * t64 + t23 * t63) * pkin(2), 0, 0, t16, -t17, 0, -t5, t4, t8, t6, t10, -t75, 0, t11 + (-t79 * qJD(5) - t5) * t68 + t70 * t65, t70 * t68 + t79 * t95 + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t37, 0.2e1 * t36, t55, t51, 0, 0, 0, 0.2e1 * t84, 0.2e1 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, 0, t59, 0, 0, 0, 0, 0, t17, t16, 0, 0, 0, 0, 0, -t75, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, 0, -t5, t4, t8, t6, t10, -t75, 0, t11 + (-pkin(4) * t16 - pkin(8) * t17) * t65 + (-t5 + (-pkin(4) * t28 + pkin(8) * t74) * qJD(5)) * t68, pkin(4) * t76 + t75 * pkin(8) + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t36, t55, t51, 0, 0, 0, t84 - t90, -t89 + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t51, 0, 0, 0, -0.2e1 * t90, -0.2e1 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t77, t17, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t95, 0, t65 * t36 - t40 * t60, t68 * t36 + t40 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t95, 0, -pkin(8) * t60, pkin(8) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
