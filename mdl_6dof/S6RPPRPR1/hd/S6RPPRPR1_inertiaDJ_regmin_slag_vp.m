% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:02
% EndTime: 2019-03-09 01:40:04
% DurationCPUTime: 0.73s
% Computational Cost: add. (1228->117), mult. (2846->216), div. (0->0), fcn. (2814->10), ass. (0->78)
t66 = sin(pkin(11));
t68 = cos(pkin(11));
t71 = sin(qJ(6));
t73 = cos(qJ(6));
t110 = -t71 * t66 + t73 * t68;
t69 = cos(pkin(10));
t105 = cos(qJ(4));
t88 = qJD(4) * t105;
t67 = sin(pkin(10));
t72 = sin(qJ(4));
t97 = t72 * t67;
t44 = qJD(4) * t97 - t69 * t88;
t93 = t66 ^ 2 + t68 ^ 2;
t112 = t93 * t44;
t111 = t110 * qJD(6);
t109 = 0.2e1 * t111;
t108 = -0.2e1 * t44;
t51 = t105 * t67 + t72 * t69;
t45 = t51 * qJD(4);
t107 = t45 * pkin(4);
t58 = sin(pkin(9)) * pkin(1) + qJ(3);
t106 = pkin(7) + t58;
t104 = t45 * t66;
t90 = t105 * t69;
t49 = -t90 + t97;
t103 = t49 * t45;
t102 = t51 * t66;
t101 = t66 * t44;
t100 = t68 * t44;
t99 = t68 * t45;
t95 = pkin(8) + qJ(5);
t46 = t106 * t67;
t47 = t106 * t69;
t21 = t46 * t88 - qJD(3) * t90 + (qJD(3) * t67 + qJD(4) * t47) * t72;
t26 = t44 * qJ(5) - t51 * qJD(5) + t107;
t7 = -t68 * t21 + t66 * t26;
t52 = -cos(pkin(9)) * pkin(1) - t69 * pkin(3) - pkin(2);
t32 = t49 * pkin(4) - t51 * qJ(5) + t52;
t37 = t105 * t47 - t72 * t46;
t14 = t66 * t32 + t68 * t37;
t6 = t66 * t21 + t68 * t26;
t13 = t68 * t32 - t66 * t37;
t87 = t93 * qJD(5);
t86 = 0.2e1 * t87;
t85 = 0.2e1 * (t67 ^ 2 + t69 ^ 2) * qJD(3);
t5 = -t68 * t51 * pkin(8) + t49 * pkin(5) + t13;
t8 = -pkin(8) * t102 + t14;
t84 = t73 * t5 - t71 * t8;
t83 = t71 * t5 + t73 * t8;
t82 = t6 * t68 + t7 * t66;
t81 = -t6 * t66 + t7 * t68;
t36 = t105 * t46 + t72 * t47;
t80 = t13 * t66 - t14 * t68;
t22 = qJD(3) * t51 + qJD(4) * t37;
t79 = t22 * t51 - t36 * t44;
t50 = t73 * t66 + t71 * t68;
t78 = -t111 * t49 - t50 * t45;
t43 = t50 * qJD(6);
t77 = -t110 * t45 + t49 * t43;
t53 = t95 * t66;
t54 = t95 * t68;
t76 = -t73 * t53 - t71 * t54;
t75 = -t71 * t53 + t73 * t54;
t74 = pkin(4) * t44 - qJ(5) * t45 - qJD(5) * t49;
t61 = -t68 * pkin(5) - pkin(4);
t35 = t110 * t51;
t34 = t50 * t51;
t28 = -qJD(5) * t50 - qJD(6) * t75;
t27 = -qJD(5) * t110 - qJD(6) * t76;
t23 = pkin(5) * t102 + t36;
t15 = -pkin(5) * t101 + t22;
t12 = t111 * t51 - t44 * t50;
t11 = t110 * t44 + t43 * t51;
t4 = pkin(8) * t101 + t7;
t3 = t45 * pkin(5) + pkin(8) * t100 + t6;
t2 = -qJD(6) * t83 + t73 * t3 - t71 * t4;
t1 = -qJD(6) * t84 - t71 * t3 - t73 * t4;
t9 = [0, 0, 0, 0, 0, 0, t85, t58 * t85, t51 * t108, 0.2e1 * t44 * t49 - 0.2e1 * t51 * t45, 0, 0, 0, 0.2e1 * t52 * t45, t52 * t108, 0.2e1 * t13 * t45 + 0.2e1 * t6 * t49 + 0.2e1 * t66 * t79, -0.2e1 * t14 * t45 - 0.2e1 * t7 * t49 + 0.2e1 * t68 * t79, -0.2e1 * t82 * t51 + 0.2e1 * (t13 * t68 + t14 * t66) * t44, 0.2e1 * t13 * t6 + 0.2e1 * t14 * t7 + 0.2e1 * t36 * t22, -0.2e1 * t35 * t11, 0.2e1 * t11 * t34 - 0.2e1 * t35 * t12, -0.2e1 * t11 * t49 + 0.2e1 * t35 * t45, -0.2e1 * t49 * t12 - 0.2e1 * t45 * t34, 0.2e1 * t103, 0.2e1 * t23 * t12 + 0.2e1 * t15 * t34 + 0.2e1 * t2 * t49 + 0.2e1 * t45 * t84, 0.2e1 * t1 * t49 - 0.2e1 * t23 * t11 + 0.2e1 * t15 * t35 - 0.2e1 * t45 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t49 + t36 * t45 + t44 * t80 + t51 * t81, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t112 * t51 + 0.2e1 * t103, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t44, t99, -t104, t112, t82, 0, 0, 0, 0, 0, -t77, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, -t22, t21, -t22 * t68 + t66 * t74, t22 * t66 + t68 * t74, t81, -t22 * pkin(4) + qJ(5) * t81 - qJD(5) * t80, -t11 * t50 + t111 * t35, -t11 * t110 - t111 * t34 - t50 * t12 - t35 * t43, -t78, -t77, 0, -t110 * t15 + t61 * t12 + t23 * t43 + t28 * t49 + t45 * t76, -t61 * t11 + t111 * t23 + t15 * t50 + t27 * t49 - t45 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t44, -t99, t104, -t112, -qJ(5) * t112 + t51 * t87 - t107, 0, 0, 0, 0, 0, t77, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, qJ(5) * t86, t50 * t109, 0.2e1 * t110 * t111 - 0.2e1 * t50 * t43, 0, 0, 0, 0.2e1 * t61 * t43, t61 * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t100, 0, t22, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, t45, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t43, 0, t28, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
