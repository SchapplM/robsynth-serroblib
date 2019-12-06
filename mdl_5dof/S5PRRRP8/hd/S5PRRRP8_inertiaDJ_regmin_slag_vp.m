% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP8
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:35
% EndTime: 2019-12-05 17:00:40
% DurationCPUTime: 0.93s
% Computational Cost: add. (637->141), mult. (1827->268), div. (0->0), fcn. (1553->8), ass. (0->92)
t42 = sin(pkin(5));
t46 = sin(qJ(2));
t104 = t42 * t46;
t43 = cos(pkin(5));
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t23 = t48 * t104 + t43 * t45;
t49 = cos(qJ(2));
t103 = t42 * t49;
t77 = qJD(2) * t103;
t13 = t23 * qJD(3) + t45 * t77;
t22 = t45 * t104 - t43 * t48;
t47 = cos(qJ(4));
t37 = qJD(4) * t47;
t44 = sin(qJ(4));
t5 = t13 * t44 + t22 * t37;
t114 = -0.4e1 * t45;
t106 = t45 * pkin(8);
t65 = -t48 * pkin(3) - t106;
t33 = -pkin(2) + t65;
t108 = pkin(7) * t48;
t36 = t47 * t108;
t113 = t44 * t33 + t36;
t40 = t47 ^ 2;
t99 = t44 ^ 2 - t40;
t69 = t99 * qJD(4);
t107 = pkin(8) * t48;
t64 = pkin(3) * t45 - t107;
t29 = t64 * qJD(3);
t92 = t45 * qJD(3);
t72 = t47 * t92;
t93 = qJD(4) * t48;
t81 = t44 * t93;
t55 = t72 + t81;
t9 = t55 * pkin(7) - t44 * t29 - t33 * t37;
t61 = pkin(4) * t44 - qJ(5) * t47;
t56 = pkin(7) + t61;
t20 = t56 * t45;
t21 = t61 * qJD(4) - t44 * qJD(5);
t62 = t47 * pkin(4) + t44 * qJ(5);
t30 = -pkin(3) - t62;
t112 = (-t30 * t48 + t106) * qJD(3) - qJD(4) * t20 - t21 * t45;
t111 = t62 * qJD(4) - t47 * qJD(5);
t76 = t44 * t92;
t10 = pkin(7) * t76 - qJD(4) * t113 + t47 * t29;
t110 = 0.2e1 * qJD(5);
t109 = pkin(7) * t44;
t39 = t45 ^ 2;
t98 = -t48 ^ 2 + t39;
t97 = qJD(2) * t46;
t96 = qJD(3) * t47;
t95 = qJD(4) * t44;
t94 = qJD(4) * t45;
t90 = t48 * qJD(3);
t89 = t48 * qJD(5);
t88 = qJ(5) * qJD(3);
t87 = -0.2e1 * pkin(2) * qJD(3);
t86 = -0.2e1 * pkin(3) * qJD(4);
t85 = t44 * t103;
t84 = pkin(4) * t92;
t83 = pkin(8) * t95;
t82 = pkin(8) * t37;
t80 = t47 * t93;
t78 = t42 * t97;
t75 = t44 * t90;
t74 = t44 * t37;
t73 = t45 * t90;
t71 = t47 * t90;
t70 = t45 * t88;
t68 = t98 * qJD(3);
t67 = 0.2e1 * t73;
t66 = t44 * t71;
t14 = t47 * t103 + t23 * t44;
t15 = t23 * t47 - t85;
t60 = t14 * t47 - t15 * t44;
t16 = -t48 * qJ(5) + t113;
t17 = -t47 * t33 + (pkin(4) + t109) * t48;
t59 = -t16 * t44 + t17 * t47;
t6 = -t13 * t47 + t22 * t95;
t12 = -qJD(3) * t22 + t48 * t77;
t2 = -qJD(4) * t85 + t12 * t44 + t23 * t37 - t47 * t78;
t54 = -t14 * t92 + t2 * t48 + t22 * t75 + t5 * t45;
t3 = -t14 * qJD(4) + t12 * t47 + t44 * t78;
t51 = t60 * qJD(4) + t2 * t44 + t3 * t47;
t4 = -t9 + t70 - t89;
t7 = -t10 - t84;
t50 = t59 * qJD(4) + t4 * t47 + t7 * t44;
t35 = pkin(8) * t80;
t24 = -t44 * t94 + t71;
t8 = t111 * t45 + t56 * t90;
t1 = (t22 * t96 + t3) * t48 + (-qJD(3) * t15 - t6) * t45;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t22 * t13 + 0.2e1 * t14 * t2 + 0.2e1 * t15 * t3; 0, 0, -t78, -t77, 0, 0, 0, 0, 0, (-t48 * t97 - t49 * t92) * t42, (t45 * t97 - t49 * t90) * t42, 0, 0, 0, 0, 0, t54, t1, t54, t60 * t90 + (t2 * t47 - t3 * t44 + (-t14 * t44 - t15 * t47) * qJD(4)) * t45, -t1, t13 * t20 + t14 * t7 + t15 * t4 + t3 * t16 + t2 * t17 + t22 * t8; 0, 0, 0, 0, t67, -0.2e1 * t68, 0, 0, 0, t45 * t87, t48 * t87, -0.2e1 * t39 * t74 + 0.2e1 * t40 * t73, t66 * t114 + 0.2e1 * t39 * t69, 0.2e1 * t45 * t81 + 0.2e1 * t98 * t96, -0.2e1 * t44 * t68 + 0.2e1 * t45 * t80, -0.2e1 * t73, 0.2e1 * t33 * t72 - 0.2e1 * t10 * t48 + 0.2e1 * (t39 * t37 + t44 * t73) * pkin(7), -0.2e1 * t9 * t48 - 0.2e1 * t113 * t92 + 0.2e1 * (-t39 * t95 + t47 * t67) * pkin(7), 0.2e1 * (qJD(3) * t20 * t44 + t7) * t48 + 0.2e1 * (-qJD(3) * t17 + t20 * t37 + t8 * t44) * t45, 0.2e1 * t59 * t90 + 0.2e1 * (-t4 * t44 + t47 * t7 + (-t16 * t47 - t17 * t44) * qJD(4)) * t45, 0.2e1 * (-t20 * t96 - t4) * t48 + 0.2e1 * (qJD(3) * t16 + t20 * t95 - t8 * t47) * t45, 0.2e1 * t16 * t4 + 0.2e1 * t17 * t7 + 0.2e1 * t20 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, 0, 0, 0, 0, 0, t6, t5, t6, t51, -t5, t51 * pkin(8) + t13 * t30 + t22 * t21; 0, 0, 0, 0, 0, 0, t90, -t92, 0, -pkin(7) * t90, pkin(7) * t92, -t45 * t69 + t66, t74 * t114 - t99 * t90, t76 - t80, t55, 0, t35 + (-pkin(3) * t47 + t109) * t94 + (t65 * t44 - t36) * qJD(3), (pkin(7) * t45 * t47 + t64 * t44) * qJD(4) + (t44 * t108 + t65 * t47) * qJD(3), t35 + (t30 * t94 - t8) * t47 - t112 * t44, t50, (-t8 + (t30 * t45 + t107) * qJD(4)) * t44 + t112 * t47, t50 * pkin(8) + t20 * t21 + t8 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t74, -0.2e1 * t69, 0, 0, 0, t44 * t86, t47 * t86, -0.2e1 * t21 * t47 + 0.2e1 * t30 * t95, 0, -0.2e1 * t21 * t44 - 0.2e1 * t30 * t37, 0.2e1 * t30 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t3, -t2, 0, t3, -t2 * pkin(4) + t3 * qJ(5) + t15 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t45 * t37 - t75, t92, t10, t9, t10 + 0.2e1 * t84, (-pkin(4) * t90 - qJ(5) * t94) * t47 + (-t48 * t88 + (pkin(4) * qJD(4) - qJD(5)) * t45) * t44, -t9 + 0.2e1 * t70 - 0.2e1 * t89, -t7 * pkin(4) + t4 * qJ(5) + t16 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t95, 0, -t82, t83, -t82, -t111, -t83, -t111 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, qJ(5) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t24, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
