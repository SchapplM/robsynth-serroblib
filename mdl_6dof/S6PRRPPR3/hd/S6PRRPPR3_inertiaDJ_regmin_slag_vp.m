% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:10:57
% EndTime: 2019-03-08 21:10:59
% DurationCPUTime: 0.70s
% Computational Cost: add. (402->127), mult. (1111->231), div. (0->0), fcn. (906->8), ass. (0->87)
t101 = pkin(8) - qJ(5);
t48 = sin(qJ(3));
t33 = t48 * qJ(4);
t51 = cos(qJ(3));
t100 = -pkin(3) * t51 - t33;
t50 = cos(qJ(6));
t42 = t50 ^ 2;
t47 = sin(qJ(6));
t94 = t47 ^ 2 - t42;
t71 = t94 * qJD(6);
t45 = cos(pkin(6));
t44 = sin(pkin(6));
t49 = sin(qJ(2));
t97 = t44 * t49;
t80 = t48 * t97;
t16 = -t45 * t51 + t80;
t17 = t45 * t48 + t51 * t97;
t52 = cos(qJ(2));
t91 = qJD(2) * t52;
t75 = t44 * t91;
t7 = -qJD(3) * t80 + (qJD(3) * t45 + t75) * t51;
t8 = qJD(3) * t17 + t48 * t75;
t56 = t8 * t48 + t7 * t51 + (t16 * t51 - t17 * t48) * qJD(3);
t84 = t48 * qJD(3);
t79 = pkin(8) * t84;
t82 = qJ(5) * qJD(3);
t13 = qJD(5) * t51 - t48 * t82 + t79;
t46 = qJ(4) + pkin(5);
t53 = -pkin(3) - pkin(4);
t39 = -pkin(9) + t53;
t98 = t39 * t48;
t99 = qJD(6) * (t46 * t51 + t98) + t13;
t54 = 2 * qJD(4);
t96 = t44 * t52;
t31 = t51 * qJD(3);
t95 = qJ(4) * t31 + qJD(4) * t48;
t43 = t51 ^ 2;
t93 = t48 ^ 2 - t43;
t92 = qJD(2) * t49;
t90 = qJD(3) * t47;
t89 = qJD(3) * t50;
t32 = qJD(6) * t47;
t88 = qJD(6) * t50;
t87 = qJD(6) * t51;
t86 = t17 * qJD(6);
t25 = t101 * t51;
t85 = t25 * qJD(6);
t83 = t51 * qJD(4);
t81 = -0.2e1 * pkin(2) * qJD(3);
t21 = -pkin(2) + t100;
t78 = t47 * t87;
t77 = t50 * t87;
t76 = t44 * t92;
t74 = t47 * t88;
t73 = t48 * t31;
t72 = t50 * t84;
t20 = pkin(4) * t51 - t21;
t70 = t93 * qJD(3);
t69 = t47 * t72;
t24 = t101 * t48;
t9 = pkin(5) * t48 + pkin(9) * t51 + t20;
t67 = t24 * t50 + t47 * t9;
t66 = t24 * t47 - t50 * t9;
t63 = qJ(4) * t7 + qJD(4) * t17;
t62 = -t16 * t47 + t50 * t96;
t61 = t16 * t50 + t47 * t96;
t60 = -t7 * t47 - t50 * t86;
t59 = -t47 * t86 + t7 * t50;
t58 = -0.2e1 * t44 ^ 2 * t49 * t91 + 0.2e1 * t16 * t8 + 0.2e1 * t17 * t7;
t57 = qJD(3) * t100 + t83;
t55 = -t83 - t85 + (-t39 * t51 + t46 * t48) * qJD(3);
t37 = qJ(4) * t54;
t29 = pkin(8) * t31;
t26 = 0.2e1 * t73;
t19 = -t31 * t47 - t48 * t88;
t18 = -t31 * t50 + t32 * t48;
t15 = pkin(3) * t84 - t95;
t14 = -qJD(5) * t48 - t51 * t82 + t29;
t12 = (t51 * t92 + t52 * t84) * t44;
t11 = -t31 * t96 + t48 * t76;
t10 = t53 * t84 + t95;
t6 = (pkin(5) * t51 + t98) * qJD(3) + t95;
t4 = qJD(6) * t62 - t47 * t76 + t8 * t50;
t3 = -qJD(6) * t61 - t8 * t47 - t50 * t76;
t2 = -qJD(6) * t67 - t47 * t14 + t50 * t6;
t1 = qJD(6) * t66 - t50 * t14 - t47 * t6;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t76, -t75, 0, 0, 0, 0, 0, -t12, t11, -t12, t56, -t11 (-t15 * t52 + t21 * t92) * t44 + t56 * pkin(8), -t11, t12, -t56, -t17 * t13 + t16 * t14 + t8 * t24 + t7 * t25 + (t10 * t52 - t20 * t92) * t44, 0, 0, 0, 0, 0 (t17 * t90 + t3) * t48 + (qJD(3) * t62 + t60) * t51 (t17 * t89 - t4) * t48 + (-qJD(3) * t61 - t59) * t51; 0, 0, 0, 0, t26, -0.2e1 * t70, 0, 0, 0, t48 * t81, t51 * t81, -0.2e1 * t15 * t51 + 0.2e1 * t21 * t84, 0, -0.2e1 * t15 * t48 - 0.2e1 * t21 * t31, 0.2e1 * t21 * t15, 0.2e1 * t10 * t48 + 0.2e1 * t20 * t31, -0.2e1 * t10 * t51 + 0.2e1 * t20 * t84, 0.2e1 * t13 * t51 - 0.2e1 * t14 * t48 + 0.2e1 * (-t24 * t51 + t25 * t48) * qJD(3), 0.2e1 * t10 * t20 - 0.2e1 * t13 * t25 + 0.2e1 * t14 * t24, -0.2e1 * t42 * t73 - 0.2e1 * t43 * t74, 0.2e1 * t43 * t71 + 0.4e1 * t51 * t69, 0.2e1 * t48 * t78 + 0.2e1 * t89 * t93, -0.2e1 * t47 * t70 + 0.2e1 * t48 * t77, t26, 0.2e1 * (t25 * t90 + t2) * t48 + 0.2e1 * (-qJD(3) * t66 + t13 * t47 - t50 * t85) * t51, 0.2e1 * (t25 * t89 + t1) * t48 + 0.2e1 * (-qJD(3) * t67 + t13 * t50 + t47 * t85) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t7, -t8, 0, t7, -pkin(3) * t8 + t63, t7, t8, 0, t53 * t8 + t63, 0, 0, 0, 0, 0, t59, t60; 0, 0, 0, 0, 0, 0, t31, -t84, 0, -t29, t79, -t29, t57, -t79, t57 * pkin(8), -t13, t14, -t83 + (-t51 * t53 + t33) * qJD(3), -qJ(4) * t13 + qJD(4) * t25 + t14 * t53, -t51 * t71 - t69, -0.4e1 * t51 * t74 + t84 * t94, t19, t18, 0, t47 * t55 - t50 * t99, t47 * t99 + t50 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t37, t54, 0, 0, t37, 0.2e1 * t74, -0.2e1 * t71, 0, 0, 0, 0.2e1 * qJD(4) * t50 - 0.2e1 * t32 * t46, -0.2e1 * qJD(4) * t47 - 0.2e1 * t46 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t29, 0, 0, -t31, t14, 0, 0, 0, 0, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t84, 0, t10, 0, 0, 0, 0, 0, -t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 + t78, -t47 * t84 + t77, t31, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t32, 0, -t39 * t88, t39 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
