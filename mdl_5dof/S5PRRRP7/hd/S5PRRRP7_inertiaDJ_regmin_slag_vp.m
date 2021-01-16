% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:23
% EndTime: 2021-01-15 16:45:30
% DurationCPUTime: 0.89s
% Computational Cost: add. (636->152), mult. (1860->292), div. (0->0), fcn. (1566->8), ass. (0->101)
t47 = sin(qJ(3));
t113 = -0.4e1 * t47;
t50 = cos(qJ(3));
t112 = t47 * qJD(5) + (pkin(7) * qJD(4) + qJ(5) * qJD(3)) * t50;
t111 = 0.2e1 * qJD(4);
t46 = sin(qJ(4));
t110 = pkin(7) * t46;
t49 = cos(qJ(4));
t109 = t49 * pkin(4);
t39 = qJD(4) * t49;
t85 = t50 * qJD(3);
t24 = t47 * t39 + t46 * t85;
t79 = pkin(7) * t85;
t17 = t24 * pkin(4) + t79;
t108 = t17 * t46;
t107 = t17 * t49;
t99 = -qJ(5) - pkin(8);
t32 = t99 * t46;
t106 = t32 * t47;
t33 = t99 * t49;
t105 = t33 * t47;
t44 = sin(pkin(5));
t48 = sin(qJ(2));
t104 = t44 * t48;
t51 = cos(qJ(2));
t103 = t44 * t51;
t102 = t46 * t50;
t101 = t47 * t49;
t100 = t49 * t50;
t62 = -t50 * pkin(3) - t47 * pkin(8);
t58 = -pkin(2) + t62;
t54 = qJD(4) * t58;
t61 = pkin(3) * t47 - pkin(8) * t50;
t55 = t61 * qJD(3);
t98 = -t46 * t55 - t49 * t54;
t87 = t47 * qJD(3);
t71 = t46 * t87;
t97 = pkin(7) * t71 + t49 * t55;
t36 = pkin(7) * t100;
t96 = t46 * t58 + t36;
t40 = t46 ^ 2;
t42 = t49 ^ 2;
t95 = t40 - t42;
t41 = t47 ^ 2;
t94 = -t50 ^ 2 + t41;
t93 = qJ(5) * t47;
t92 = qJD(2) * t48;
t91 = qJD(3) * t46;
t90 = qJD(3) * t49;
t89 = qJD(4) * t46;
t88 = qJD(4) * t50;
t84 = -0.2e1 * pkin(2) * qJD(3);
t83 = -0.2e1 * pkin(3) * qJD(4);
t82 = t49 * t103;
t81 = pkin(4) * t87;
t80 = pkin(4) * t89;
t78 = t47 * t89;
t77 = t46 * t88;
t76 = t49 * t88;
t45 = cos(pkin(5));
t21 = t47 * t104 - t45 * t50;
t75 = t21 * t89;
t30 = (pkin(4) * t46 + pkin(7)) * t47;
t74 = t30 * t89;
t73 = t44 * t92;
t72 = qJD(2) * t103;
t70 = t46 * t39;
t69 = t47 * t85;
t68 = t49 * t85;
t37 = -pkin(3) - t109;
t67 = -t37 + t109;
t66 = t95 * qJD(4);
t65 = t94 * qJD(3);
t64 = 0.2e1 * t69;
t63 = t46 * t68;
t60 = pkin(4) * t40 + t37 * t49;
t22 = t50 * t104 + t45 * t47;
t14 = -t22 * t46 - t82;
t56 = t46 * t103 - t22 * t49;
t59 = -t14 * t49 + t46 * t56;
t13 = t22 * qJD(3) + t47 * t72;
t7 = t13 * t46 + t21 * t39;
t8 = -t13 * t49 + t75;
t23 = t68 - t78;
t53 = t49 * t87 + t77;
t52 = qJ(5) * t78 - t112 * t49 - t46 * t54 + t97;
t29 = t49 * t58;
t20 = -t46 * qJD(5) + t99 * t39;
t19 = -t49 * qJD(5) - t99 * t89;
t16 = -t46 * t93 + t96;
t12 = -qJD(3) * t21 + t50 * t72;
t11 = -t49 * t93 + t29 + (-pkin(4) - t110) * t50;
t10 = -t96 * qJD(4) + t97;
t9 = t53 * pkin(7) + t98;
t6 = (pkin(7) * qJD(3) + qJ(5) * qJD(4)) * t101 + t112 * t46 + t98;
t5 = -qJD(4) * t82 + t12 * t49 - t22 * t89 + t46 * t73;
t4 = t56 * qJD(4) - t12 * t46 + t49 * t73;
t3 = t52 + t81;
t2 = (t21 * t90 + t5) * t50 + (qJD(3) * t56 - t8) * t47;
t1 = (t21 * t91 - t4) * t50 + (qJD(3) * t14 + t7) * t47;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t21 * t13 + 0.2e1 * t14 * t4 - 0.2e1 * t5 * t56; 0, 0, -t73, -t72, 0, 0, 0, 0, 0, (-t50 * t92 - t51 * t87) * t44, (t47 * t92 - t51 * t85) * t44, 0, 0, 0, 0, 0, t1, t2, t1, t2, t59 * t85 + (-t4 * t49 - t46 * t5 + (t14 * t46 + t49 * t56) * qJD(4)) * t47, t4 * t11 + t13 * t30 + t14 * t3 + t5 * t16 + t21 * t17 + t56 * t6; 0, 0, 0, 0, t64, -0.2e1 * t65, 0, 0, 0, t47 * t84, t50 * t84, -0.2e1 * t41 * t70 + 0.2e1 * t42 * t69, t95 * t41 * t111 + t113 * t63, 0.2e1 * t47 * t77 + 0.2e1 * t94 * t90, -0.2e1 * t46 * t65 + 0.2e1 * t47 * t76, -0.2e1 * t69, 0.2e1 * t29 * t87 - 0.2e1 * t10 * t50 + 0.2e1 * (t41 * t39 + t46 * t69) * pkin(7), -0.2e1 * t9 * t50 - 0.2e1 * t96 * t87 + 0.2e1 * (-t41 * t89 + t49 * t64) * pkin(7), 0.2e1 * (t30 * t91 - t3) * t50 + 0.2e1 * (qJD(3) * t11 + t30 * t39 + t108) * t47, 0.2e1 * (t30 * t90 - t6) * t50 + 0.2e1 * (-qJD(3) * t16 + t107 - t74) * t47, 0.2e1 * (-t11 * t49 - t16 * t46) * t85 + 0.2e1 * (-t3 * t49 + t46 * t6 + (t11 * t46 - t16 * t49) * qJD(4)) * t47, 0.2e1 * t11 * t3 - 0.2e1 * t16 * t6 + 0.2e1 * t30 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, 0, 0, 0, 0, 0, t8, t7, t8, t7, t59 * qJD(4) - t4 * t46 + t5 * t49, pkin(4) * t75 + t13 * t37 + t14 * t20 + t19 * t56 + t4 * t32 - t5 * t33; 0, 0, 0, 0, 0, 0, t85, -t87, 0, -t79, pkin(7) * t87, -t47 * t66 + t63, t113 * t70 - t95 * t85, t71 - t76, t53, 0, (pkin(8) * t100 + (-pkin(3) * t49 + t110) * t47) * qJD(4) + (t62 * t46 - t36) * qJD(3), (pkin(7) * t101 + t61 * t46) * qJD(4) + (pkin(7) * t102 + t62 * t49) * qJD(3), -t107 - t20 * t50 + (t37 * t102 + t106) * qJD(3) + (t30 * t46 + t60 * t47) * qJD(4), t108 - t19 * t50 + (t37 * t100 + t105) * qJD(3) + (t67 * t47 * t46 + t30 * t49) * qJD(4), (-t32 * t85 - t20 * t47 - t6 + (-t11 + t105) * qJD(4)) * t49 + (t33 * t85 + t19 * t47 - t3 + (-t16 + t106) * qJD(4)) * t46, pkin(4) * t74 + t11 * t20 - t16 * t19 + t17 * t37 + t3 * t32 + t6 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t70, -0.2e1 * t66, 0, 0, 0, t46 * t83, t49 * t83, -0.2e1 * t67 * t89, t60 * t111, -0.2e1 * t19 * t49 - 0.2e1 * t20 * t46 + 0.2e1 * (-t32 * t49 + t33 * t46) * qJD(4), 0.2e1 * t33 * t19 + 0.2e1 * t32 * t20 + 0.2e1 * t37 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t5, t4, -t5, 0, t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t24, t87, t10, t9, t52 + 0.2e1 * t81, t6, -t23 * pkin(4), t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t89, 0, -pkin(8) * t39, pkin(8) * t89, t20, t19, -pkin(4) * t39, t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t23, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t39, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
