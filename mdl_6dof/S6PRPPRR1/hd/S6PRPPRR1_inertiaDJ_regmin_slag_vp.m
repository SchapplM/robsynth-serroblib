% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:14
% EndTime: 2019-03-08 19:16:16
% DurationCPUTime: 0.53s
% Computational Cost: add. (587->103), mult. (1654->197), div. (0->0), fcn. (1716->12), ass. (0->85)
t56 = cos(qJ(6));
t46 = t56 ^ 2;
t53 = sin(qJ(6));
t83 = t53 ^ 2 - t46;
t71 = t83 * qJD(6);
t47 = sin(pkin(12));
t50 = cos(pkin(12));
t54 = sin(qJ(5));
t95 = cos(qJ(5));
t37 = t95 * t47 + t54 * t50;
t32 = t37 * qJD(5);
t97 = 0.2e1 * t32;
t48 = sin(pkin(11));
t41 = t48 * pkin(2) + qJ(4);
t96 = pkin(8) + t41;
t33 = t96 * t47;
t34 = t96 * t50;
t17 = -t54 * t33 + t95 * t34;
t10 = qJD(4) * t37 + qJD(5) * t17;
t94 = t10 * t53;
t93 = t10 * t56;
t49 = sin(pkin(6));
t51 = cos(pkin(11));
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t30 = (t48 * t57 + t51 * t55) * t49;
t25 = qJD(2) * t30;
t61 = t48 * t55 - t51 * t57;
t29 = t61 * t49;
t92 = t29 * t25;
t72 = qJD(5) * t95;
t88 = t54 * t47;
t31 = qJD(5) * t88 - t50 * t72;
t91 = t37 * t31;
t90 = t53 * t31;
t89 = t53 * t32;
t87 = t56 * t31;
t86 = t56 * t32;
t74 = t95 * t50;
t36 = -t74 + t88;
t85 = -t36 * t87 + t37 * t86;
t84 = t47 ^ 2 + t50 ^ 2;
t82 = qJD(2) * t49;
t81 = qJD(6) * t53;
t80 = qJD(6) * t56;
t79 = -0.2e1 * pkin(5) * qJD(6);
t78 = t53 * t87;
t77 = t37 * t81;
t76 = t53 * t80;
t75 = -t51 * pkin(2) - pkin(3);
t26 = t61 * t82;
t73 = t84 * t26;
t70 = 0.2e1 * t84 * qJD(4);
t69 = pkin(5) * t31 - pkin(9) * t32;
t68 = pkin(5) * t37 + pkin(9) * t36;
t52 = cos(pkin(6));
t22 = -t30 * t47 + t52 * t50;
t23 = t30 * t50 + t52 * t47;
t8 = t54 * t22 + t95 * t23;
t67 = t29 * t56 - t53 * t8;
t66 = t29 * t53 + t56 * t8;
t38 = -t50 * pkin(4) + t75;
t11 = t36 * pkin(5) - t37 * pkin(9) + t38;
t65 = t56 * t11 - t53 * t17;
t64 = t53 * t11 + t56 * t17;
t63 = -t22 * t47 + t23 * t50;
t62 = t31 * t36 - t37 * t32;
t6 = qJD(5) * t8 - t26 * t37;
t7 = -t22 * t95 + t54 * t23;
t60 = t6 * t53 + t7 * t80;
t59 = -t6 * t56 + t7 * t81;
t58 = t37 * t80 - t90;
t13 = t77 + t87;
t14 = t36 * t80 + t89;
t35 = t37 ^ 2;
t18 = t32 * pkin(5) + t31 * pkin(9);
t16 = t33 * t95 + t54 * t34;
t12 = t36 * t81 - t86;
t9 = t33 * t72 - qJD(4) * t74 + (qJD(4) * t47 + qJD(5) * t34) * t54;
t5 = -t22 * t72 + t26 * t74 + (qJD(5) * t23 - t26 * t47) * t54;
t4 = -qJD(6) * t64 + t56 * t18 + t53 * t9;
t3 = -qJD(6) * t65 - t53 * t18 + t56 * t9;
t2 = -qJD(6) * t66 + t25 * t56 + t53 * t5;
t1 = -qJD(6) * t67 - t25 * t53 + t56 * t5;
t15 = [0, 0, 0, 0, -0.2e1 * t30 * t26 + 0.2e1 * t92, 0, 0, 0, -0.2e1 * t26 * t63 + 0.2e1 * t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t55 * t82, -t57 * t82 (-t25 * t51 - t26 * t48) * pkin(2), -t25 * t50, t25 * t47, -t73, qJD(4) * t63 + t25 * t75 - t41 * t73, 0, 0, 0, 0, 0, t25 * t36 + t29 * t32, t25 * t37 - t29 * t31, 0, 0, 0, 0, 0, t2 * t36 + t32 * t67 + t37 * t60 - t7 * t90, t1 * t36 - t32 * t66 - t37 * t59 - t7 * t87; 0, 0, 0, 0, 0, 0, 0, t70, t41 * t70, -0.2e1 * t91, 0.2e1 * t62, 0, 0, 0, t38 * t97, -0.2e1 * t38 * t31, -0.2e1 * t35 * t76 - 0.2e1 * t46 * t91, 0.2e1 * t35 * t71 + 0.4e1 * t37 * t78, -0.2e1 * t36 * t77 + 0.2e1 * t85, -0.2e1 * t36 * t58 - 0.2e1 * t37 * t89, t36 * t97, 0.2e1 * t16 * t58 + 0.2e1 * t32 * t65 + 0.2e1 * t4 * t36 + 0.2e1 * t37 * t94, -0.2e1 * t13 * t16 + 0.2e1 * t3 * t36 - 0.2e1 * t32 * t64 + 0.2e1 * t37 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t62 + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, 0, 0, 0, 0, -t12, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, t59, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t32, 0, -t10, t9, -t37 * t71 - t78, t31 * t83 - 0.4e1 * t37 * t76, t14, -t12, 0, -t93 + t69 * t53 + (t16 * t53 - t56 * t68) * qJD(6), t94 + t69 * t56 + (t16 * t56 + t53 * t68) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t31, 0, 0, 0, 0, 0, t12, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t76, -0.2e1 * t71, 0, 0, 0, t53 * t79, t56 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t58, t32, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t81, 0, -pkin(9) * t80, pkin(9) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t15;
