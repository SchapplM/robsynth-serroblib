% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x21]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:29
% EndTime: 2019-03-08 19:58:32
% DurationCPUTime: 0.88s
% Computational Cost: add. (859->158), mult. (2362->286), div. (0->0), fcn. (2218->10), ass. (0->100)
t53 = sin(qJ(4));
t106 = -0.4e1 * t53;
t52 = sin(qJ(5));
t56 = cos(qJ(4));
t86 = t56 * qJD(4);
t55 = cos(qJ(5));
t89 = qJD(5) * t55;
t60 = t52 * t86 + t53 * t89;
t105 = t60 * pkin(5);
t50 = cos(pkin(11));
t41 = -t50 * pkin(2) - pkin(3);
t69 = -t56 * pkin(4) - t53 * pkin(9);
t32 = t41 + t69;
t21 = t52 * t32;
t101 = t55 * t56;
t48 = sin(pkin(11));
t40 = t48 * pkin(2) + pkin(8);
t34 = t40 * t101;
t104 = -t21 - t34;
t44 = t52 ^ 2;
t46 = t55 ^ 2;
t72 = (t44 - t46) * qJD(5);
t103 = t40 * t52;
t102 = t53 * t55;
t100 = -qJ(6) - pkin(9);
t68 = pkin(4) * t53 - pkin(9) * t56;
t35 = t68 * qJD(4);
t99 = -t32 * t89 - t52 * t35;
t43 = t53 * qJD(4);
t76 = t40 * t43;
t98 = t55 * t35 + t52 * t76;
t45 = t53 ^ 2;
t96 = -t56 ^ 2 + t45;
t95 = qJ(6) * t53;
t94 = qJ(6) * t55;
t49 = sin(pkin(6));
t93 = qJD(2) * t49;
t54 = sin(qJ(2));
t57 = cos(qJ(2));
t19 = (t48 * t57 + t50 * t54) * t49;
t51 = cos(pkin(6));
t13 = t19 * t53 - t51 * t56;
t92 = qJD(4) * t13;
t91 = qJD(4) * t55;
t90 = qJD(5) * t52;
t88 = qJD(5) * t56;
t87 = t55 * qJD(6);
t85 = -0.2e1 * pkin(4) * qJD(5);
t84 = 0.2e1 * qJD(4) * t41;
t83 = pkin(5) * t90;
t82 = t52 * t88;
t81 = t55 * t88;
t80 = t13 * t90;
t79 = t44 * t86;
t78 = t52 * t89;
t77 = t53 * t86;
t75 = t55 * t86;
t74 = t40 * t86;
t73 = qJD(5) * t100;
t71 = t96 * qJD(4);
t70 = t52 * t75;
t14 = t19 * t56 + t51 * t53;
t63 = t48 * t54 - t50 * t57;
t18 = t63 * t49;
t10 = t14 * t55 + t18 * t52;
t9 = -t14 * t52 + t18 * t55;
t67 = t10 * t55 - t52 * t9;
t66 = -t10 * t52 - t55 * t9;
t22 = t55 * t32;
t11 = -t53 * t94 + t22 + (-pkin(5) - t103) * t56;
t12 = -t52 * t95 - t104;
t65 = -t11 * t55 - t12 * t52;
t64 = t11 * t52 - t12 * t55;
t17 = t63 * t93;
t7 = t14 * qJD(4) - t17 * t53;
t62 = t13 * t89 + t7 * t52;
t61 = -t7 * t55 + t80;
t26 = t53 * t90 - t75;
t27 = t55 * t43 + t82;
t16 = qJD(2) * t19;
t8 = -t17 * t56 - t92;
t1 = -t10 * qJD(5) + t16 * t55 - t8 * t52;
t2 = t9 * qJD(5) + t16 * t52 + t8 * t55;
t59 = t66 * qJD(5) - t1 * t52 + t2 * t55;
t23 = t52 * t73 + t87;
t24 = -t52 * qJD(6) + t55 * t73;
t37 = t100 * t52;
t38 = t100 * t55;
t58 = t23 * t55 - t24 * t52 + (-t37 * t55 + t38 * t52) * qJD(5);
t42 = -t55 * pkin(5) - pkin(4);
t39 = t46 * t86;
t36 = t46 * t77;
t29 = t52 * t43 - t81;
t25 = (pkin(5) * t52 + t40) * t53;
t15 = t74 + t105;
t6 = qJD(5) * t104 + t98;
t5 = t27 * t40 + t99;
t4 = (-qJ(6) * qJD(5) - qJD(4) * t40) * t102 + (-qJD(6) * t53 + (-qJ(6) * qJD(4) - qJD(5) * t40) * t56) * t52 - t99;
t3 = -t53 * t87 + (pkin(5) * t53 - t56 * t94) * qJD(4) + (-t34 + (-t32 + t95) * t52) * qJD(5) + t98;
t20 = [0, 0, 0, 0, 0.2e1 * t18 * t16 - 0.2e1 * t19 * t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t13 * t7; 0, 0, -t54 * t93, -t57 * t93 (-t16 * t50 - t17 * t48) * pkin(2), 0, 0, 0, 0, 0, -t16 * t56 + t18 * t43, t16 * t53 + t18 * t86, 0, 0, 0, 0, 0 (t52 * t92 - t1) * t56 + (qJD(4) * t9 + t62) * t53 (t13 * t91 + t2) * t56 + (-qJD(4) * t10 - t61) * t53, t66 * t86 + (-t67 * qJD(5) - t1 * t55 - t2 * t52) * t53, t1 * t11 + t10 * t4 + t2 * t12 + t13 * t15 + t7 * t25 + t9 * t3; 0, 0, 0, 0, 0, 0.2e1 * t77, -0.2e1 * t71, 0, 0, 0, t53 * t84, t56 * t84, -0.2e1 * t45 * t78 + 0.2e1 * t36, t106 * t70 + 0.2e1 * t45 * t72, 0.2e1 * t53 * t82 + 0.2e1 * t96 * t91, -0.2e1 * t52 * t71 + 0.2e1 * t53 * t81, -0.2e1 * t77, 0.2e1 * t22 * t43 - 0.2e1 * t6 * t56 + 0.2e1 * (t45 * t89 + t52 * t77) * t40, -0.2e1 * t45 * t40 * t90 - 0.2e1 * t5 * t56 + 0.2e1 * (-t21 + t34) * t43, 0.2e1 * t65 * t86 + 0.2e1 * (t64 * qJD(5) - t3 * t55 - t4 * t52) * t53, 0.2e1 * t11 * t3 + 0.2e1 * t12 * t4 + 0.2e1 * t25 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t67 * qJD(4) - t7) * t56 + (t59 + t92) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t64 * qJD(4) - t15) * t56 + (qJD(4) * t25 + t65 * qJD(5) - t3 * t52 + t4 * t55) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t36 + 0.2e1 * (t44 - 0.1e1) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, 0, 0, 0, 0, t61, t62, t59, pkin(5) * t80 + t1 * t37 + t10 * t23 - t2 * t38 + t9 * t24 + t7 * t42; 0, 0, 0, 0, 0, 0, 0, t86, -t43, 0, -t74, t76, -t53 * t72 + t70, t106 * t78 + t39 - t79, t29, t27, 0 (pkin(9) * t101 + (-pkin(4) * t55 + t103) * t53) * qJD(5) + (t69 * t52 - t34) * qJD(4) (t40 * t102 + t68 * t52) * qJD(5) + (t56 * t103 + t69 * t55) * qJD(4) (-t37 * t86 - t24 * t53 + t4 + (t38 * t53 - t11) * qJD(5)) * t55 + (t38 * t86 - t23 * t53 - t3 + (t37 * t53 - t12) * qJD(5)) * t52, t11 * t24 + t12 * t23 + t15 * t42 + t25 * t83 + t3 * t37 - t4 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t86, 0, 0, 0, 0, 0, -t27, t29, t39 + t79 (-t83 + (-t37 * t52 - t38 * t55) * qJD(4)) * t56 + (qJD(4) * t42 + t58) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t78, -0.2e1 * t72, 0, 0, 0, t52 * t85, t55 * t85, 0.2e1 * t58, -0.2e1 * t38 * t23 + 0.2e1 * t37 * t24 + 0.2e1 * t42 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t60, t43, t6, t5, t26 * pkin(5), t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t26, 0, -t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -t90, 0, -pkin(9) * t89, pkin(9) * t90, -pkin(5) * t89, t24 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t20;
