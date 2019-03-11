% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x20]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:38
% EndTime: 2019-03-08 18:43:40
% DurationCPUTime: 0.64s
% Computational Cost: add. (588->112), mult. (1990->235), div. (0->0), fcn. (2139->14), ass. (0->96)
t58 = cos(qJ(6));
t45 = t58 ^ 2;
t55 = sin(qJ(6));
t103 = t55 ^ 2 - t45;
t76 = t103 * qJD(6);
t52 = cos(pkin(12));
t53 = cos(pkin(7));
t105 = t52 * t53;
t49 = sin(pkin(7));
t57 = sin(qJ(3));
t107 = t49 * t57;
t48 = sin(pkin(12));
t50 = sin(pkin(6));
t54 = cos(pkin(6));
t60 = cos(qJ(3));
t21 = t50 * (t105 * t57 + t48 * t60) + t54 * t107;
t47 = sin(pkin(13));
t41 = pkin(3) * t47 + pkin(9);
t108 = t41 * t55;
t51 = cos(pkin(13));
t106 = t51 * t60;
t59 = cos(qJ(5));
t104 = t58 * t59;
t56 = sin(qJ(5));
t44 = t56 ^ 2;
t102 = -t59 ^ 2 + t44;
t101 = qJD(3) * t49;
t100 = qJD(3) * t57;
t99 = qJD(5) * t55;
t98 = qJD(5) * t58;
t97 = qJD(6) * t55;
t96 = qJD(6) * t58;
t95 = qJD(6) * t59;
t29 = (t47 * t60 + t51 * t57) * t49;
t22 = t29 * t56 - t53 * t59;
t94 = t22 * qJD(6);
t93 = t56 * qJD(5);
t92 = t59 * qJD(5);
t91 = -0.2e1 * pkin(5) * qJD(6);
t90 = t60 * t105;
t88 = t59 * t108;
t87 = t41 * t104;
t42 = -pkin(3) * t51 - pkin(4);
t86 = 0.2e1 * qJD(5) * t42;
t85 = t55 * t95;
t84 = t58 * t95;
t83 = t60 * t101;
t82 = t55 * t93;
t81 = t55 * t96;
t80 = t56 * t92;
t79 = t41 * t93;
t78 = t58 * t93;
t77 = t58 * t92;
t75 = t102 * qJD(5);
t74 = t56 * t77;
t73 = -pkin(5) * t59 - pkin(10) * t56;
t72 = pkin(5) * t56 - pkin(10) * t59;
t62 = t54 * t49 * t60 + (-t48 * t57 + t90) * t50;
t12 = t21 * t47 - t51 * t62;
t13 = t51 * t21 + t47 * t62;
t30 = -t49 * t50 * t52 + t53 * t54;
t8 = t13 * t59 + t30 * t56;
t71 = t12 * t58 - t55 * t8;
t70 = t12 * t55 + t58 * t8;
t23 = t29 * t59 + t53 * t56;
t28 = -t106 * t49 + t107 * t47;
t69 = t23 * t58 + t28 * t55;
t68 = t23 * t55 - t28 * t58;
t19 = -t54 * t83 + (-qJD(3) * t90 + t100 * t48) * t50;
t61 = t21 * qJD(3);
t11 = -t51 * t19 - t47 * t61;
t4 = qJD(5) * t8 + t56 * t11;
t7 = t13 * t56 - t30 * t59;
t67 = t4 * t55 + t7 * t96;
t66 = -t4 * t58 + t7 * t97;
t26 = (-t47 * t57 + t106) * t101;
t17 = qJD(5) * t23 + t56 * t26;
t64 = t17 * t55 + t58 * t94;
t63 = -t17 * t58 + t55 * t94;
t37 = t72 * qJD(5);
t35 = t42 + t73;
t34 = t82 - t84;
t33 = -t55 * t92 - t56 * t96;
t32 = t78 + t85;
t31 = t56 * t97 - t77;
t25 = qJD(3) * t29;
t16 = -t26 * t59 + t29 * t93 - t53 * t92;
t15 = t55 * t79 + t58 * t37 + (-t35 * t55 - t87) * qJD(6);
t14 = t41 * t78 - t55 * t37 + (-t35 * t58 + t88) * qJD(6);
t10 = -t19 * t47 + t51 * t61;
t6 = -qJD(6) * t69 + t55 * t16 + t58 * t25;
t5 = qJD(6) * t68 + t58 * t16 - t55 * t25;
t3 = -t11 * t59 + t13 * t93 - t30 * t92;
t2 = -qJD(6) * t70 + t58 * t10 + t55 * t3;
t1 = -qJD(6) * t71 - t55 * t10 + t58 * t3;
t9 = [0, 0, 0, 0, 0, 0.2e1 * t10 * t12 + 0.2e1 * t11 * t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t10 * t28 + t11 * t29 + t12 * t25 + t13 * t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0.2e1 * t25 * t28 + 0.2e1 * t26 * t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t61, t19 (-t10 * t51 + t11 * t47) * pkin(3), 0, 0, 0, 0, 0, -t10 * t59 + t12 * t93, t10 * t56 + t12 * t92, 0, 0, 0, 0, 0 (t7 * t99 - t2) * t59 + (qJD(5) * t71 + t67) * t56 (t7 * t98 - t1) * t59 + (-qJD(5) * t70 - t66) * t56; 0, 0, 0, -t49 * t100, -t83 (-t25 * t51 + t26 * t47) * pkin(3), 0, 0, 0, 0, 0, -t25 * t59 + t28 * t93, t25 * t56 + t28 * t92, 0, 0, 0, 0, 0 (t22 * t99 - t6) * t59 + (-qJD(5) * t68 + t64) * t56 (t22 * t98 - t5) * t59 + (-qJD(5) * t69 - t63) * t56; 0, 0, 0, 0, 0, 0, 0.2e1 * t80, -0.2e1 * t75, 0, 0, 0, t56 * t86, t59 * t86, -0.2e1 * t44 * t81 + 0.2e1 * t45 * t80, 0.2e1 * t44 * t76 - 0.4e1 * t55 * t74, 0.2e1 * t102 * t98 + 0.2e1 * t56 * t85, -0.2e1 * t55 * t75 + 0.2e1 * t56 * t84, -0.2e1 * t80, 0.2e1 * t35 * t78 - 0.2e1 * t15 * t59 + 0.2e1 * (t44 * t96 + t55 * t80) * t41, -0.2e1 * t35 * t82 - 0.2e1 * t14 * t59 + 0.2e1 * (-t44 * t97 + t74) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, 0, 0, 0, 0, 0, t66, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, 0, 0, 0, 0, t63, t64; 0, 0, 0, 0, 0, 0, 0, 0, t92, -t93, 0, -t41 * t92, t79, t55 * t77 - t56 * t76, -t103 * t92 - 0.4e1 * t56 * t81, t34, t32, 0 (pkin(10) * t104 + (-pkin(5) * t58 + t108) * t56) * qJD(6) + (t55 * t73 - t87) * qJD(5) (t41 * t56 * t58 + t55 * t72) * qJD(6) + (t58 * t73 + t88) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t92, 0, 0, 0, 0, 0, -t32, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t81, -0.2e1 * t76, 0, 0, 0, t55 * t91, t58 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t33, t93, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t97, 0, -pkin(10) * t96, pkin(10) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
