% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:03
% EndTime: 2019-03-08 19:45:05
% DurationCPUTime: 0.69s
% Computational Cost: add. (730->133), mult. (1965->248), div. (0->0), fcn. (1969->10), ass. (0->89)
t102 = cos(qJ(4));
t47 = sin(pkin(11));
t49 = cos(pkin(11));
t52 = sin(qJ(4));
t33 = t102 * t47 + t52 * t49;
t28 = t33 * qJD(4);
t79 = t102 * t49;
t95 = t52 * t47;
t32 = -t79 + t95;
t48 = sin(pkin(6));
t55 = cos(qJ(2));
t53 = sin(qJ(2));
t88 = qJD(2) * t53;
t107 = (-t28 * t55 + t32 * t88) * t48;
t77 = qJD(4) * t102;
t27 = qJD(4) * t95 - t49 * t77;
t106 = (t27 * t55 + t33 * t88) * t48;
t51 = sin(qJ(6));
t45 = t51 ^ 2;
t54 = cos(qJ(6));
t90 = -t54 ^ 2 + t45;
t76 = t90 * qJD(6);
t93 = pkin(8) + qJ(3);
t35 = t93 * t47;
t36 = t93 * t49;
t12 = (qJD(3) * t47 + qJD(4) * t36) * t52 - qJD(3) * t79 + t35 * t77;
t105 = -0.2e1 * t27;
t104 = 2 * qJD(5);
t103 = pkin(4) + pkin(9);
t101 = t32 * t51;
t100 = t32 * t54;
t99 = t48 ^ 2 * t53;
t98 = t48 * t53;
t97 = t48 * t55;
t96 = t51 * t28;
t94 = t54 * t28;
t91 = t47 ^ 2 + t49 ^ 2;
t89 = qJD(2) * t48;
t87 = qJD(2) * t55;
t86 = qJD(6) * t51;
t85 = qJD(6) * t54;
t84 = qJD(6) * t103;
t83 = qJ(5) * qJD(6);
t82 = t51 * t94;
t37 = t48 * t88;
t81 = t48 * t87;
t80 = t51 * t85;
t41 = -t49 * pkin(3) - pkin(2);
t78 = t91 * t55;
t75 = 0.2e1 * t91 * qJD(3);
t22 = t102 * t35 + t52 * t36;
t68 = -t33 * qJ(5) + t41;
t11 = t103 * t32 + t68;
t14 = t33 * pkin(5) + t22;
t74 = t54 * t11 + t51 * t14;
t73 = t51 * t11 - t54 * t14;
t50 = cos(pkin(6));
t25 = -t47 * t98 + t50 * t49;
t26 = t50 * t47 + t49 * t98;
t72 = -t25 * t47 + t26 * t49;
t71 = t27 * qJ(5) - t33 * qJD(5);
t70 = -qJ(5) * t28 - qJD(5) * t32;
t16 = -t102 * t25 + t52 * t26;
t67 = -t51 * t16 + t54 * t97;
t66 = t54 * t16 + t51 * t97;
t63 = -t51 * t27 + t33 * t85;
t62 = t54 * t27 + t33 * t86;
t61 = t32 * t85 + t96;
t60 = t32 * t86 - t94;
t17 = t102 * t26 + t52 * t25;
t59 = t102 * t36 - t52 * t35;
t8 = -t28 * pkin(5) - t12;
t58 = t8 + (qJ(5) * t32 + t103 * t33) * qJD(6);
t15 = -t32 * pkin(5) + t59;
t57 = -qJD(6) * t15 - t103 * t27 - t70;
t13 = qJD(3) * t33 + qJD(4) * t59;
t30 = t32 ^ 2;
t21 = t33 * t105;
t20 = t32 * pkin(4) + t68;
t10 = t28 * pkin(4) + t71;
t9 = -t27 * pkin(5) + t13;
t7 = t103 * t28 + t71;
t6 = qJD(4) * t17 + t33 * t81;
t5 = -t25 * t77 - t79 * t81 + (qJD(4) * t26 + t47 * t81) * t52;
t4 = qJD(6) * t66 + t54 * t37 + t51 * t6;
t3 = t67 * qJD(6) - t51 * t37 + t54 * t6;
t2 = -t74 * qJD(6) - t51 * t7 + t54 * t9;
t1 = t73 * qJD(6) - t51 * t9 - t54 * t7;
t18 = [0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t72 * t48 - t99) * t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t16 * t6 - 0.2e1 * t17 * t5 - 0.2e1 * t87 * t99, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t37, -t81, -t49 * t37, t47 * t37, t78 * t89, t72 * qJD(3) + (-pkin(2) * t53 + qJ(3) * t78) * t89, 0, 0, 0, 0, 0, t107, t106, -t16 * t27 - t17 * t28 + t5 * t32 + t6 * t33, -t107, -t106, -t17 * t12 + t16 * t13 + t6 * t22 - t5 * t59 + (-t10 * t55 + t20 * t88) * t48, 0, 0, 0, 0, 0, t5 * t100 + t60 * t17 - t27 * t66 + t3 * t33, -t5 * t101 + t61 * t17 - t67 * t27 - t4 * t33; 0, 0, 0, 0, 0, 0, t75, qJ(3) * t75, t21, 0.2e1 * t27 * t32 - 0.2e1 * t33 * t28, 0, 0, 0, 0.2e1 * t41 * t28, t41 * t105, 0.2e1 * t12 * t32 + 0.2e1 * t13 * t33 - 0.2e1 * t22 * t27 - 0.2e1 * t28 * t59, -0.2e1 * t10 * t32 - 0.2e1 * t20 * t28, -0.2e1 * t10 * t33 + 0.2e1 * t20 * t27, 0.2e1 * t20 * t10 - 0.2e1 * t12 * t59 + 0.2e1 * t22 * t13, 0.2e1 * t45 * t32 * t28 + 0.2e1 * t30 * t80, -0.2e1 * t30 * t76 + 0.4e1 * t32 * t82, 0.2e1 * t32 * t63 + 0.2e1 * t33 * t96, -0.2e1 * t32 * t62 + 0.2e1 * t33 * t94, t21, -0.2e1 * t8 * t100 + 0.2e1 * t60 * t15 + 0.2e1 * t2 * t33 + 0.2e1 * t73 * t27, 0.2e1 * t1 * t33 + 0.2e1 * t8 * t101 + 0.2e1 * t61 * t15 + 0.2e1 * t74 * t27; 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, 0, -t28, t27, t10, 0, 0, 0, 0, 0, -t63, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, t6, -t5, -t6 * pkin(4) - t5 * qJ(5) + t17 * qJD(5), 0, 0, 0, 0, 0, t17 * t85 - t5 * t51, -t17 * t86 - t5 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t28, 0, -t13, t12, pkin(4) * t27 + t70, t13, -t12, -t13 * pkin(4) - t12 * qJ(5) + qJD(5) * t59, -t32 * t76 + t82, -t90 * t28 - 0.4e1 * t32 * t80, -t62, -t63, 0, t51 * t58 - t54 * t57, t51 * t57 + t54 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, qJ(5) * t104, -0.2e1 * t80, 0.2e1 * t76, 0, 0, 0, 0.2e1 * qJD(5) * t51 + 0.2e1 * t54 * t83, 0.2e1 * qJD(5) * t54 - 0.2e1 * t51 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, 0, t13, 0, 0, 0, 0, 0, -t62, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t60, -t27, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85, 0, t51 * t84, t54 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t18;
