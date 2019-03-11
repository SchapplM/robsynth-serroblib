% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:21
% EndTime: 2019-03-09 01:49:23
% DurationCPUTime: 0.87s
% Computational Cost: add. (559->119), mult. (1299->225), div. (0->0), fcn. (1095->6), ass. (0->83)
t52 = sin(qJ(4));
t47 = sin(pkin(9));
t48 = cos(pkin(9));
t51 = sin(qJ(6));
t53 = cos(qJ(6));
t100 = -t51 * t47 + t53 * t48;
t97 = t100 * qJD(6);
t102 = t52 * t97;
t54 = cos(qJ(4));
t101 = t97 * t54;
t49 = -pkin(7) + qJ(2);
t83 = t54 * qJD(4);
t75 = t49 * t83;
t85 = t52 * qJD(2);
t99 = t75 + t85;
t29 = t53 * t47 + t51 * t48;
t98 = t29 * t54;
t96 = 0.2e1 * t97;
t95 = t48 * pkin(8);
t94 = t52 * pkin(4);
t24 = t29 * qJD(6);
t93 = t24 * t52;
t92 = t49 * t52;
t50 = pkin(1) + qJ(3);
t89 = pkin(8) + qJ(5);
t68 = -t54 * qJ(5) + t94;
t31 = t68 + t50;
t38 = t48 * t92;
t16 = t47 * t31 + t38;
t88 = t47 ^ 2 + t48 ^ 2;
t87 = qJD(5) * t52;
t46 = t54 ^ 2;
t86 = t46 * qJD(2);
t43 = t52 * qJD(4);
t84 = t54 * qJD(2);
t82 = qJ(2) * qJD(2);
t81 = t47 * t92;
t22 = -t54 * qJD(5) + qJD(3) + (pkin(4) * t54 + qJ(5) * t52) * qJD(4);
t12 = t47 * t22 + t99 * t48;
t78 = t47 * t43;
t77 = t48 * t43;
t76 = t52 * t83;
t74 = -t47 * t49 + pkin(5);
t73 = pkin(5) * t47 - t49;
t72 = t88 * t54;
t71 = t88 * qJD(5);
t70 = 0.2e1 * t76;
t69 = 0.2e1 * t71;
t18 = t48 * t22;
t11 = -t47 * t99 + t18;
t67 = -t11 * t48 - t12 * t47;
t66 = -t11 * t47 + t12 * t48;
t26 = t48 * t31;
t13 = t74 * t52 - t54 * t95 + t26;
t14 = -t47 * t54 * pkin(8) + t16;
t65 = t53 * t13 - t51 * t14;
t64 = t51 * t13 + t53 * t14;
t15 = t26 - t81;
t63 = -t15 * t47 + t16 * t48;
t36 = t89 * t47;
t37 = t89 * t48;
t62 = -t53 * t36 - t51 * t37;
t61 = -t51 * t36 + t53 * t37;
t60 = t29 * t83 + t102;
t59 = -t100 * t83 + t93;
t21 = t100 * t54;
t30 = t49 * t43 - t84;
t56 = qJD(4) * t100;
t55 = 0.2e1 * qJD(2);
t42 = -t48 * pkin(5) - pkin(4);
t27 = t73 * t54;
t19 = -t73 * t43 - t84;
t10 = -t51 * t77 - t53 * t78 + t101;
t9 = -qJD(4) * t98 - t102;
t8 = -qJD(6) * t98 - t52 * t56;
t7 = -t54 * t56 + t93;
t6 = -qJD(5) * t29 - qJD(6) * t61;
t5 = -qJD(5) * t100 - qJD(6) * t62;
t4 = pkin(8) * t78 + t12;
t3 = -t47 * t85 + t18 + (t52 * t95 + t74 * t54) * qJD(4);
t2 = -t64 * qJD(6) + t53 * t3 - t51 * t4;
t1 = -t65 * qJD(6) - t51 * t3 - t53 * t4;
t17 = [0, 0, 0, 0, t55, 0.2e1 * t82, t55, 0.2e1 * qJD(3), 0.2e1 * t50 * qJD(3) + 0.2e1 * t82, -0.2e1 * t76, 0.2e1 * (t52 ^ 2 - t46) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t52 + 0.2e1 * t50 * t83, 0.2e1 * qJD(3) * t54 - 0.2e1 * t50 * t43, -0.2e1 * t47 * t86 + 0.2e1 * t11 * t52 + 0.2e1 * (t15 + 0.2e1 * t81) * t83, -0.2e1 * t48 * t86 - 0.2e1 * t12 * t52 + 0.2e1 * (-t16 + 0.2e1 * t38) * t83, 0.2e1 * t67 * t54 + 0.2e1 * (t15 * t48 + t16 * t47) * t43, 0.2e1 * t15 * t11 + 0.2e1 * t16 * t12 + 0.2e1 * (-t52 * t75 + t86) * t49, 0.2e1 * t21 * t8, -0.2e1 * t21 * t10 - 0.2e1 * t8 * t98, 0.2e1 * t21 * t83 + 0.2e1 * t8 * t52, -0.2e1 * t10 * t52 - 0.2e1 * t83 * t98, t70, 0.2e1 * t27 * t10 + 0.2e1 * t19 * t98 + 0.2e1 * t2 * t52 + 0.2e1 * t65 * t83, 0.2e1 * t1 * t52 + 0.2e1 * t19 * t21 + 0.2e1 * t27 * t8 - 0.2e1 * t64 * t83; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, 0, 0, 0, -t83, t43, -t48 * t83, t47 * t83, -t88 * t43, t67, 0, 0, 0, 0, 0, t59, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86 + t66 * t52 + (t63 - 0.2e1 * t92) * t83, 0, 0, 0, 0, 0, -t54 * t10 + t9 * t52, t7 * t52 - t54 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t88) * t70, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t83, 0, -t30, -t99, t48 * t84 - t47 * t87 + (t68 * t47 - t38) * qJD(4), -t47 * t84 - t48 * t87 + (t68 * t48 + t81) * qJD(4), t66, -t30 * pkin(4) + t66 * qJ(5) + t63 * qJD(5), t21 * t97 + t8 * t29, -t29 * t10 + t100 * t8 - t21 * t24 - t97 * t98, t60, -t59, 0, t42 * t10 - t100 * t19 + t27 * t24 + t6 * t52 + t62 * t83, t19 * t29 + t27 * t97 + t42 * t8 + t5 * t52 - t61 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t83, -t77, t78, qJD(4) * t72, t52 * t71 + (qJ(5) * t72 - t94) * qJD(4), 0, 0, 0, 0, 0, -t100 * t43 - t54 * t24, t29 * t43 - t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, qJ(5) * t69, t29 * t96, 0.2e1 * t100 * t97 - 0.2e1 * t29 * t24, 0, 0, 0, 0.2e1 * t42 * t24, t42 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t77, 0, t30, 0, 0, 0, 0, 0, t10, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t10, t83, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, -t24, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t17;
