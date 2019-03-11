% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:32
% EndTime: 2019-03-09 03:32:35
% DurationCPUTime: 0.96s
% Computational Cost: add. (838->148), mult. (1586->242), div. (0->0), fcn. (1113->4), ass. (0->94)
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t99 = t52 * qJ(4);
t65 = t50 * pkin(3) - t99;
t32 = qJ(2) + t65;
t27 = t50 * pkin(8) + t32;
t54 = -pkin(1) - pkin(7);
t34 = (pkin(4) - t54) * t52;
t49 = sin(qJ(5));
t51 = cos(qJ(5));
t115 = t51 * t27 + t49 * t34;
t46 = t50 ^ 2;
t48 = t52 ^ 2;
t72 = (t46 - t48) * qJD(3);
t45 = t49 ^ 2;
t47 = t51 ^ 2;
t102 = t45 - t47;
t71 = t102 * qJD(5);
t53 = -pkin(3) - pkin(8);
t104 = t52 * t53;
t63 = t49 * pkin(5) - t51 * qJ(6);
t31 = qJ(4) + t63;
t112 = t63 * qJD(5) - t49 * qJD(6);
t88 = t52 * qJD(3);
t39 = t54 * t88;
t64 = pkin(5) * t51 + qJ(6) * t49;
t57 = -pkin(4) - t64;
t6 = t112 * t50 + t57 * t88 + t39;
t114 = (t31 * t50 - t104) * qJD(5) + t6;
t105 = t50 * t53;
t90 = t50 * qJD(4);
t41 = t50 * t54;
t33 = -t50 * pkin(4) + t41;
t92 = t33 * qJD(5);
t113 = (t99 + t105) * qJD(3) + t90 - t92;
t16 = t65 * qJD(3) - t90;
t43 = t50 * qJD(3);
t76 = pkin(3) * t88 + qJ(4) * t43 + qJD(2);
t11 = (qJD(3) * pkin(8) - qJD(4)) * t52 + t76;
t38 = t54 * t43;
t29 = -pkin(4) * t43 + t38;
t5 = -qJD(5) * t115 - t49 * t11 + t51 * t29;
t111 = -t64 * qJD(5) + t51 * qJD(6);
t95 = qJD(5) * t51;
t96 = qJD(5) * t49;
t4 = -t51 * t11 + t27 * t96 - t49 * t29 - t34 * t95;
t75 = qJ(6) * t43;
t87 = t52 * qJD(6);
t2 = -t4 - t75 + t87;
t84 = pkin(5) * t43;
t3 = -t5 + t84;
t7 = t52 * qJ(6) + t115;
t61 = -t49 * t27 + t51 * t34;
t8 = -t52 * pkin(5) - t61;
t67 = t49 * t8 + t51 * t7;
t1 = t67 * qJD(5) + t2 * t49 - t3 * t51;
t110 = 0.2e1 * qJD(2);
t109 = 0.2e1 * qJD(4);
t108 = 0.2e1 * qJD(6);
t107 = t31 * t52;
t12 = qJD(4) - t111;
t106 = t50 * t12;
t101 = t45 + t47;
t98 = qJD(3) * t49;
t97 = qJD(3) * t51;
t94 = qJD(5) * t52;
t93 = qJD(5) * t53;
t86 = qJ(2) * qJD(3);
t85 = qJ(4) * qJD(5);
t83 = t49 * t94;
t82 = t49 * t93;
t81 = t51 * t94;
t80 = t51 * t93;
t79 = t51 * t88;
t78 = t49 * t95;
t77 = t50 * t88;
t74 = t101 * t50;
t70 = qJD(5) * (t46 + t48);
t69 = t49 * t79;
t66 = t49 * t7 - t51 * t8;
t30 = -pkin(4) * t88 + t39;
t56 = t30 + (qJ(4) * t50 - t104) * qJD(5);
t10 = t57 * t50 + t41;
t55 = qJD(5) * t10 - t106 + (-t105 - t107) * qJD(3);
t35 = -0.2e1 * t77;
t25 = -t49 * t43 + t81;
t24 = t51 * t70;
t23 = t49 * t88 + t50 * t95;
t22 = t51 * t43 + t83;
t21 = t49 * t70;
t20 = t50 * t96 - t79;
t19 = qJD(3) * t74;
t13 = -t52 * qJD(4) + t76;
t9 = [0, 0, 0, 0, t110, qJ(2) * t110, t35, 0.2e1 * t72, 0, 0, 0, 0.2e1 * qJD(2) * t50 + 0.2e1 * t52 * t86, 0.2e1 * qJD(2) * t52 - 0.2e1 * t50 * t86, 0, -0.2e1 * t13 * t50 - 0.2e1 * t32 * t88, -0.2e1 * t13 * t52 + 0.2e1 * t32 * t43, 0.2e1 * t32 * t13, 0.2e1 * t45 * t77 + 0.2e1 * t46 * t78, -0.2e1 * t46 * t71 + 0.4e1 * t50 * t69, -0.2e1 * t49 * t72 + 0.2e1 * t50 * t81, -0.2e1 * t50 * t83 - 0.2e1 * t51 * t72, t35, 0.2e1 * (-t33 * t97 + t5) * t52 + 0.2e1 * (-t61 * qJD(3) - t30 * t51 + t49 * t92) * t50, 0.2e1 * (t33 * t98 + t4) * t52 + 0.2e1 * (qJD(3) * t115 + t30 * t49 + t51 * t92) * t50, 0.2e1 * (-t10 * t97 - t3) * t52 + 0.2e1 * (qJD(3) * t8 + t10 * t96 - t6 * t51) * t50, 0.2e1 * t67 * t88 + 0.2e1 * (-t66 * qJD(5) + t2 * t51 + t3 * t49) * t50, 0.2e1 * (-t10 * t98 + t2) * t52 + 0.2e1 * (-qJD(3) * t7 - t10 * t95 - t6 * t49) * t50, 0.2e1 * t10 * t6 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t24, t21, 0, -t24 (t66 * qJD(3) + t6) * t50 + (qJD(3) * t10 - t1) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (0.1e1 - t101) * t77; 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t88, 0, -t38, -t39, t16, t38, t39, -t16 * t54, -t50 * t71 + t69, -t102 * t88 - 0.4e1 * t50 * t78, -t22, -t25, 0, -t113 * t51 + t56 * t49, t113 * t49 + t56 * t51, t114 * t49 + t55 * t51, -t1, -t114 * t51 + t55 * t49, t1 * t53 + t10 * t12 + t6 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t88, 0, t43, t88, -t16, 0, 0, 0, 0, 0, t23, -t20, t23, -t19, t20, t106 + (t53 * t74 + t107) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, qJ(4) * t109, -0.2e1 * t78, 0.2e1 * t71, 0, 0, 0, 0.2e1 * qJD(4) * t49 + 0.2e1 * t51 * t85, 0.2e1 * qJD(4) * t51 - 0.2e1 * t49 * t85, 0.2e1 * t12 * t49 + 0.2e1 * t31 * t95, 0, -0.2e1 * t12 * t51 + 0.2e1 * t31 * t96, 0.2e1 * t31 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, 0, t38, 0, 0, 0, 0, 0, -t22, -t25, -t22, 0, t25, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t20, -t43, t5, t4, t5 - 0.2e1 * t84, t111 * t50 - t63 * t88, -t4 - 0.2e1 * t75 + 0.2e1 * t87, -t3 * pkin(5) + t2 * qJ(6) + t7 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t25, t22, 0, -t25 (-qJ(6) * t94 + t84) * t51 + (t75 + (pkin(5) * qJD(5) - qJD(6)) * t52) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t95, 0, -t82, -t80, -t82, t112, t80, -t112 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t95, -t96, 0, t95, -t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, qJ(6) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t23, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
