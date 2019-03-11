% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPPRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:02
% EndTime: 2019-03-09 01:32:05
% DurationCPUTime: 1.11s
% Computational Cost: add. (1499->121), mult. (2787->206), div. (0->0), fcn. (2659->8), ass. (0->86)
t45 = sin(pkin(10));
t46 = cos(pkin(10));
t97 = sin(qJ(5));
t74 = qJD(5) * t97;
t49 = cos(qJ(5));
t87 = qJD(5) * t49;
t30 = -t45 * t74 + t46 * t87;
t32 = t45 * t49 + t46 * t97;
t22 = t32 * t30;
t106 = 0.2e1 * t22;
t29 = t32 * qJD(5);
t33 = -t45 * t97 + t46 * t49;
t23 = t33 * t29;
t114 = -0.2e1 * t23;
t113 = t106 + t114;
t38 = -cos(pkin(9)) * pkin(1) - pkin(2) - qJ(4);
t98 = -pkin(7) + t38;
t28 = t98 * t45;
t73 = t97 * qJD(4);
t76 = t49 * t98;
t83 = t49 * qJD(4);
t109 = (qJD(5) * t76 - t73) * t46 - t28 * t74 - t45 * t83;
t21 = t32 * t29;
t112 = -t33 * t30 + t21;
t39 = sin(pkin(9)) * pkin(1) + qJ(3);
t34 = pkin(4) * t45 + t39;
t53 = pkin(5) * t32 - pkin(8) * t33 + t34;
t110 = -qJD(6) * t53 - t109;
t47 = sin(qJ(6));
t43 = t47 ^ 2;
t48 = cos(qJ(6));
t44 = t48 ^ 2;
t72 = qJD(6) * (t43 - t44);
t35 = (t45 ^ 2 + t46 ^ 2) * qJD(4);
t70 = t98 * t97;
t11 = t49 * t28 + t46 * t70;
t101 = t30 * pkin(5);
t102 = t29 * pkin(8);
t61 = qJD(3) + t101 + t102;
t86 = qJD(6) * t47;
t1 = t11 * t86 + t110 * t48 - t47 * t61;
t85 = qJD(6) * t48;
t2 = -t11 * t85 + t110 * t47 + t48 * t61;
t3 = -t11 * t47 + t48 * t53;
t4 = t48 * t11 + t47 * t53;
t65 = t3 * t47 - t4 * t48;
t108 = qJD(6) * t65 + t1 * t47 - t2 * t48;
t31 = t33 ^ 2;
t105 = 0.2e1 * qJD(3);
t10 = t28 * t97 - t46 * t76;
t7 = t28 * t87 - t45 * t73 + (qJD(5) * t70 + t83) * t46;
t104 = t10 * t7;
t103 = t29 * pkin(5);
t100 = t7 * t47;
t99 = t7 * t48;
t96 = t29 * t48;
t95 = t43 * t29;
t26 = t44 * t29;
t94 = t47 * t30;
t93 = t48 * t30;
t92 = -t21 * t48 + t33 * t93;
t88 = t43 + t44;
t84 = t39 * qJD(3);
t81 = -0.2e1 * pkin(5) * qJD(6);
t80 = t47 * t96;
t79 = t33 * t86;
t78 = t47 * t85;
t77 = t32 ^ 2 + t31;
t75 = t88 * t30;
t71 = t31 * t78;
t69 = -pkin(8) * t30 + t103;
t68 = pkin(5) * t33 + pkin(8) * t32;
t66 = t3 * t48 + t4 * t47;
t64 = t10 * t29 - t7 * t33;
t63 = t10 * t30 + t7 * t32;
t60 = -t29 * t47 + t33 * t85;
t14 = t79 + t96;
t15 = t32 * t85 + t94;
t54 = -qJD(6) * t66 - t1 * t48 - t2 * t47;
t50 = t109 * t32 + t11 * t30 + t64;
t18 = t33 * t26;
t17 = t33 * t95;
t13 = t32 * t86 - t93;
t12 = t26 + t95;
t6 = t33 * t72 + t80;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, 0.2e1 * t84, 0, 0, 0, 0, 0, 0, t45 * t105, t46 * t105, 0.2e1 * t35, -0.2e1 * t35 * t38 + 0.2e1 * t84, t114, 0.2e1 * t112, 0, t106, 0, 0, 0.2e1 * qJD(3) * t32 + 0.2e1 * t30 * t34, 0.2e1 * qJD(3) * t33 - 0.2e1 * t29 * t34, -0.2e1 * t50, 0.2e1 * t34 * qJD(3) + 0.2e1 * t109 * t11 + 0.2e1 * t104, -0.2e1 * t18 - 0.2e1 * t71, 0.2e1 * t31 * t72 + 0.4e1 * t33 * t80, -0.2e1 * t32 * t79 + 0.2e1 * t92, -0.2e1 * t17 + 0.2e1 * t71, -0.2e1 * t32 * t60 - 0.2e1 * t33 * t94, t106, 0.2e1 * t10 * t60 + 0.2e1 * t100 * t33 + 0.2e1 * t2 * t32 + 0.2e1 * t3 * t30, 0.2e1 * t1 * t32 - 0.2e1 * t10 * t14 - 0.2e1 * t4 * t30 + 0.2e1 * t33 * t99, 0.2e1 * t108 * t33 + 0.2e1 * t29 * t66, -0.2e1 * t1 * t4 + 0.2e1 * t2 * t3 + 0.2e1 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t33 - t11 * t29 + t63, 0, 0, 0, 0, 0, 0, 0, t112 * t48 + t92, 0, t29 * t65 + t33 * t54 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t18 - 0.2e1 * t17 + 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t50, 0, 0, 0, 0, 0, 0, -t113 * t47 - t77 * t85, -t113 * t48 + t77 * t86, 0, -t30 * t65 + t32 * t54 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112 * (-0.1e1 + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t32 * t75 - 0.2e1 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, 0, t30, -t29, 0, qJD(3), 0, 0, 0, 0, 0, 0, -t13, -t15, t12, -t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, -t30, 0, -t7, -t109, 0, 0, -t6, -0.4e1 * t33 * t78 - t26 + t95, t15, t6, -t13, 0, -t99 + t69 * t47 + (t10 * t47 - t48 * t68) * qJD(6), t100 + t69 * t48 + (t10 * t48 + t47 * t68) * qJD(6), t54, -t7 * pkin(5) + pkin(8) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, -t12, -t102 * t88 - t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t30, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t60, t75, pkin(8) * t75 - t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t78, -0.2e1 * t72, 0, -0.2e1 * t78, 0, 0, t47 * t81, t48 * t81, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, -t60, t30, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, -t86, 0, -pkin(8) * t85, pkin(8) * t86, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
