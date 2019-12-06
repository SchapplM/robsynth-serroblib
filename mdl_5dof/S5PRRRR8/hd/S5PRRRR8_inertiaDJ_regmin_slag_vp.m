% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:41
% EndTime: 2019-12-05 17:16:44
% DurationCPUTime: 0.66s
% Computational Cost: add. (709->120), mult. (1929->229), div. (0->0), fcn. (1831->10), ass. (0->92)
t60 = cos(qJ(5));
t53 = t60 ^ 2;
t56 = sin(qJ(5));
t95 = t56 ^ 2 - t53;
t77 = t95 * qJD(5);
t107 = qJD(3) + qJD(4);
t106 = pkin(7) + pkin(8);
t57 = sin(qJ(4));
t58 = sin(qJ(3));
t61 = cos(qJ(4));
t62 = cos(qJ(3));
t39 = t57 * t58 - t61 * t62;
t25 = t107 * t39;
t40 = t57 * t62 + t61 * t58;
t105 = t40 * t25;
t104 = t40 * t56;
t103 = t40 * t60;
t54 = sin(pkin(5));
t59 = sin(qJ(2));
t102 = t54 * t59;
t63 = cos(qJ(2));
t101 = t54 * t63;
t26 = t107 * t40;
t100 = t56 * t26;
t99 = t60 * t25;
t98 = t60 * t26;
t43 = t106 * t58;
t44 = t106 * t62;
t29 = -t57 * t43 + t61 * t44;
t79 = qJD(3) * t106;
t75 = t62 * t79;
t76 = t58 * t79;
t14 = t29 * qJD(4) - t57 * t76 + t61 * t75;
t28 = t61 * t43 + t57 * t44;
t51 = qJD(5) * t60;
t97 = t14 * t56 + t28 * t51;
t49 = -t61 * pkin(3) - pkin(4);
t93 = qJD(4) * t57;
t84 = pkin(3) * t93;
t96 = t49 * t51 + t56 * t84;
t94 = qJD(2) * t59;
t92 = qJD(4) * t61;
t91 = qJD(5) * t56;
t90 = t58 * qJD(3);
t89 = t62 * qJD(3);
t88 = -0.2e1 * pkin(2) * qJD(3);
t87 = pkin(4) * t91;
t86 = pkin(4) * t51;
t85 = pkin(3) * t90;
t83 = pkin(3) * t92;
t82 = t54 * t94;
t81 = qJD(2) * t101;
t80 = t56 * t51;
t50 = -t62 * pkin(3) - pkin(2);
t78 = -0.4e1 * t56 * t103;
t19 = t39 * pkin(4) - t40 * pkin(9) + t50;
t74 = t60 * t19 - t56 * t29;
t73 = t56 * t19 + t60 * t29;
t55 = cos(pkin(5));
t31 = -t58 * t102 + t55 * t62;
t32 = t62 * t102 + t55 * t58;
t18 = t57 * t31 + t61 * t32;
t48 = t57 * pkin(3) + pkin(9);
t72 = t39 * t48 - t40 * t49;
t71 = t49 * t91 - t60 * t84;
t70 = t60 * t101 + t56 * t18;
t69 = t56 * t101 - t60 * t18;
t68 = -t56 * t25 + t40 * t51;
t67 = t40 * t91 + t99;
t66 = t39 * t91 - t98;
t65 = qJD(3) * t32 + t58 * t81;
t64 = -t25 * t49 - t26 * t48 + (-t39 * t61 + t40 * t57) * qJD(4) * pkin(3);
t46 = 0.2e1 * t80;
t38 = -0.2e1 * t77;
t37 = t40 ^ 2;
t27 = qJD(3) * t31 + t62 * t81;
t23 = t28 * t91;
t17 = -t61 * t31 + t57 * t32;
t16 = t39 * t51 + t100;
t13 = t43 * t92 + t44 * t93 + t57 * t75 + t61 * t76;
t11 = t26 * pkin(4) + t25 * pkin(9) + t85;
t10 = -t40 * t77 - t56 * t99;
t9 = qJD(5) * t78 + t95 * t25;
t8 = t18 * qJD(4) + t57 * t27 + t61 * t65;
t7 = -t61 * t27 - t31 * t92 + t32 * t93 + t57 * t65;
t6 = t17 * t91 - t8 * t60;
t5 = t17 * t51 + t8 * t56;
t4 = qJD(5) * t69 + t56 * t7 + t60 * t82;
t3 = qJD(5) * t70 - t56 * t82 + t60 * t7;
t2 = -t73 * qJD(5) + t60 * t11 + t56 * t13;
t1 = -t74 * qJD(5) - t56 * t11 + t60 * t13;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t82, -t81, 0, 0, 0, 0, 0, (-t62 * t94 - t63 * t90) * t54, (t58 * t94 - t63 * t89) * t54, 0, 0, 0, 0, 0, (-t26 * t63 + t39 * t94) * t54, (t25 * t63 + t40 * t94) * t54, 0, 0, 0, 0, 0, t8 * t104 + t17 * t68 - t26 * t70 + t4 * t39, t8 * t103 - t17 * t67 + t26 * t69 + t3 * t39; 0, 0, 0, 0, 0.2e1 * t58 * t89, 0.2e1 * (-t58 ^ 2 + t62 ^ 2) * qJD(3), 0, 0, 0, t58 * t88, t62 * t88, -0.2e1 * t105, 0.2e1 * t25 * t39 - 0.2e1 * t40 * t26, 0, 0, 0, 0.2e1 * t50 * t26 + 0.2e1 * t39 * t85, -0.2e1 * t50 * t25 + 0.2e1 * t40 * t85, -0.2e1 * t53 * t105 - 0.2e1 * t37 * t80, -t25 * t78 + 0.2e1 * t37 * t77, -0.2e1 * t39 * t67 + 0.2e1 * t40 * t98, -0.2e1 * t40 * t100 - 0.2e1 * t39 * t68, 0.2e1 * t39 * t26, 0.2e1 * t14 * t104 + 0.2e1 * t2 * t39 + 0.2e1 * t74 * t26 + 0.2e1 * t68 * t28, 0.2e1 * t1 * t39 + 0.2e1 * t14 * t103 - 0.2e1 * t73 * t26 - 0.2e1 * t67 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t27, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, t89, -t90, 0, -pkin(7) * t89, pkin(7) * t90, 0, 0, -t25, -t26, 0, -t14, t13, t10, t9, t16, -t66, 0, t23 + (-t72 * qJD(5) - t14) * t60 + t64 * t56, t64 * t60 + t72 * t91 + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t84, -0.2e1 * t83, t46, t38, 0, 0, 0, 0.2e1 * t71, 0.2e1 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, -t14, t13, t10, t9, t16, -t66, 0, t23 + (pkin(4) * t25 - pkin(9) * t26) * t56 + (-t14 + (-pkin(4) * t40 - pkin(9) * t39) * qJD(5)) * t60, pkin(4) * t67 + pkin(9) * t66 + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t83, t46, t38, 0, 0, 0, t71 - t87, -t86 + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t38, 0, 0, 0, -0.2e1 * t87, -0.2e1 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t68, t26, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t91, 0, -t48 * t51 - t56 * t83, t48 * t91 - t60 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t91, 0, -pkin(9) * t51, pkin(9) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
