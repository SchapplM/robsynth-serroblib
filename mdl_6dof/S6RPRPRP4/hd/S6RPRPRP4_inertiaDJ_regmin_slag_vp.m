% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:01
% EndTime: 2019-03-09 03:13:04
% DurationCPUTime: 0.92s
% Computational Cost: add. (853->146), mult. (1705->242), div. (0->0), fcn. (1245->6), ass. (0->96)
t38 = sin(pkin(9)) * pkin(1) + pkin(7);
t108 = pkin(4) + t38;
t52 = cos(qJ(3));
t53 = -pkin(3) - pkin(8);
t103 = t52 * t53;
t39 = -cos(pkin(9)) * pkin(1) - pkin(2);
t50 = sin(qJ(3));
t97 = t50 * qJ(4);
t57 = t39 - t97;
t15 = t57 + t103;
t29 = t108 * t50;
t49 = sin(qJ(5));
t51 = cos(qJ(5));
t113 = t15 * t51 + t29 * t49;
t47 = t52 ^ 2;
t70 = qJD(3) * (t50 ^ 2 - t47);
t44 = t49 ^ 2;
t46 = t51 ^ 2;
t101 = t44 - t46;
t69 = t101 * qJD(5);
t41 = t50 * qJD(3);
t68 = pkin(3) * t41 - t50 * qJD(4);
t98 = qJ(4) * t52;
t11 = (pkin(8) * t50 - t98) * qJD(3) + t68;
t42 = t52 * qJD(3);
t31 = t38 * t42;
t20 = pkin(4) * t42 + t31;
t5 = -qJD(5) * t113 - t49 * t11 + t20 * t51;
t93 = qJD(5) * t51;
t94 = qJD(5) * t49;
t4 = -t11 * t51 + t15 * t94 - t20 * t49 - t29 * t93;
t86 = qJ(6) * qJD(3);
t73 = t52 * t86;
t89 = t50 * qJD(6);
t2 = -t4 + t73 + t89;
t83 = pkin(5) * t42;
t3 = -t5 - t83;
t7 = qJ(6) * t50 + t113;
t60 = -t15 * t49 + t29 * t51;
t8 = -pkin(5) * t50 - t60;
t65 = t49 * t8 + t51 * t7;
t1 = qJD(5) * t65 + t2 * t49 - t3 * t51;
t112 = 0.2e1 * qJD(4);
t111 = 0.2e1 * qJD(6);
t110 = t52 * pkin(3);
t62 = -pkin(5) * t49 + qJ(6) * t51;
t21 = qJD(5) * t62 + t49 * qJD(6);
t63 = pkin(5) * t51 + qJ(6) * t49;
t6 = t21 * t52 + (-t63 - t108) * t41;
t109 = t6 * t49;
t17 = qJD(5) * t63 - t51 * qJD(6) + qJD(4);
t107 = t17 * t52;
t19 = t108 * t41;
t106 = t19 * t49;
t32 = qJ(4) - t62;
t105 = t32 * t52;
t104 = t50 * t53;
t30 = t108 * t52;
t100 = t44 + t46;
t96 = qJD(3) * t49;
t95 = qJD(3) * t51;
t92 = qJD(5) * t52;
t91 = qJD(5) * t53;
t90 = t30 * qJD(5);
t88 = t52 * qJD(4);
t87 = qJ(4) * qJD(5);
t85 = qJD(3) * qJ(4);
t84 = 0.2e1 * qJD(3) * t39;
t82 = t49 * t92;
t81 = t49 * t91;
t80 = t51 * t92;
t79 = t51 * t91;
t78 = t50 * t42;
t77 = t38 * t41;
t76 = t51 * t41;
t75 = t51 * t42;
t74 = t49 * t93;
t72 = t100 * t50;
t35 = 0.2e1 * t78;
t67 = t49 * t76;
t64 = t49 * t7 - t51 * t8;
t59 = t104 + t105;
t58 = -t98 - t104;
t56 = pkin(5) * t41 - qJ(6) * t92;
t55 = t50 * t86 + (pkin(5) * qJD(5) - qJD(6)) * t52;
t54 = t88 + (-t97 - t110) * qJD(3);
t33 = t53 * t75;
t28 = t57 - t110;
t27 = -t41 * t49 + t80;
t26 = t42 * t49 + t50 * t93;
t25 = t76 + t82;
t24 = t50 * t94 - t75;
t23 = qJD(3) * t72;
t22 = t52 * t85 - t68;
t9 = t52 * t63 + t30;
t10 = [0, 0, 0, 0, t35, -0.2e1 * t70, 0, 0, 0, t50 * t84, t52 * t84, 0, -0.2e1 * t22 * t52 - 0.2e1 * t28 * t41, 0.2e1 * t22 * t50 - 0.2e1 * t28 * t42, -0.2e1 * t28 * t22, -0.2e1 * t44 * t78 + 0.2e1 * t47 * t74, -0.2e1 * t47 * t69 - 0.4e1 * t52 * t67, 0.2e1 * t49 * t70 - 0.2e1 * t50 * t80, 0.2e1 * t50 * t82 + 0.2e1 * t51 * t70, t35, 0.2e1 * (-t30 * t95 + t5) * t50 + 0.2e1 * (qJD(3) * t60 - t19 * t51 - t49 * t90) * t52, 0.2e1 * (t30 * t96 + t4) * t50 + 0.2e1 * (-qJD(3) * t113 - t51 * t90 + t106) * t52, 0.2e1 * (-t9 * t95 - t3) * t50 + 0.2e1 * (-qJD(3) * t8 + t6 * t51 - t9 * t94) * t52, 0.2e1 * t65 * t41 + 0.2e1 * (qJD(5) * t64 - t2 * t51 - t3 * t49) * t52, 0.2e1 * (-t9 * t96 + t2) * t50 + 0.2e1 * (qJD(3) * t7 + t9 * t93 + t109) * t52, 0.2e1 * t2 * t7 + 0.2e1 * t3 * t8 + 0.2e1 * t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (qJD(3) * t64 + t6) * t50 + (qJD(3) * t9 - t1) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (0.1e1 - t100) * t35; 0, 0, 0, 0, 0, 0, t42, -t41, 0, -t31, t77, t54, t31, -t77, t54 * t38, t52 * t69 + t67, -t101 * t41 + 0.4e1 * t52 * t74, -t24, -t26, 0, -t106 + t33 + (-t50 * t85 + t88) * t51 + (t30 * t51 + t49 * t58) * qJD(5) (qJD(5) * t58 - t19) * t51 + (-t88 - t90 + (t97 - t103) * qJD(3)) * t49, t109 + t33 + (-t32 * t41 + t107) * t51 + (-t49 * t59 + t51 * t9) * qJD(5), -t1 (qJD(5) * t59 - t6) * t51 + (qJD(5) * t9 + t107 + (-t32 * t50 + t103) * qJD(3)) * t49, t1 * t53 + t9 * t17 + t6 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t42, 0, t41, t42, t22, 0, 0, 0, 0, 0, t26, -t24, t26, -t23, t24, t50 * t17 + (t53 * t72 + t105) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, qJ(4) * t112, -0.2e1 * t74, 0.2e1 * t69, 0, 0, 0, 0.2e1 * qJD(4) * t49 + 0.2e1 * t51 * t87, 0.2e1 * qJD(4) * t51 - 0.2e1 * t49 * t87, 0.2e1 * t17 * t49 + 0.2e1 * t32 * t93, 0, -0.2e1 * t17 * t51 + 0.2e1 * t32 * t94, 0.2e1 * t32 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, t31, 0, 0, 0, 0, 0, -t24, -t26, -t24, 0, t26, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t25, t42, t5, t4, t5 + 0.2e1 * t83, -t49 * t56 + t51 * t55, -t4 + 0.2e1 * t73 + 0.2e1 * t89, -pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t27, t25, 0, -t27, t49 * t55 + t51 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t93, 0, -t81, -t79, -t81, -t21, t79, t21 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t93, -t94, 0, t93, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, qJ(6) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t27, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, 0, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
