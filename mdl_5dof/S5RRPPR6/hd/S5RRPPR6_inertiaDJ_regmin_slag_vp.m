% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:05
% EndTime: 2019-12-31 19:33:08
% DurationCPUTime: 0.69s
% Computational Cost: add. (973->105), mult. (2348->223), div. (0->0), fcn. (2178->8), ass. (0->75)
t58 = sin(pkin(9));
t60 = cos(pkin(9));
t62 = sin(qJ(5));
t64 = cos(qJ(5));
t95 = -t62 * t58 + t64 * t60;
t94 = t95 * qJD(5);
t44 = t64 * t58 + t62 * t60;
t38 = t44 * qJD(5);
t93 = 0.2e1 * t94;
t59 = sin(pkin(8));
t51 = t59 * pkin(2) + qJ(4);
t92 = pkin(7) + t51;
t63 = sin(qJ(2));
t65 = cos(qJ(2));
t85 = -qJ(3) - pkin(6);
t77 = qJD(2) * t85;
t34 = t65 * qJD(3) + t63 * t77;
t61 = cos(pkin(8));
t66 = -t63 * qJD(3) + t65 * t77;
t22 = t59 * t34 - t61 * t66;
t47 = t85 * t63;
t48 = t85 * t65;
t30 = -t61 * t47 - t59 * t48;
t91 = t30 * t22;
t43 = t59 * t65 + t61 * t63;
t90 = t43 * t58;
t82 = t65 * qJD(2);
t83 = t63 * qJD(2);
t36 = -t59 * t83 + t61 * t82;
t89 = t58 * t36;
t88 = t60 * t36;
t35 = t43 * qJD(2);
t55 = pkin(2) * t83;
t17 = t35 * pkin(3) - t36 * qJ(4) - t43 * qJD(4) + t55;
t23 = t61 * t34 + t59 * t66;
t6 = t58 * t17 + t60 * t23;
t41 = t59 * t63 - t61 * t65;
t78 = -t65 * pkin(2) - pkin(1);
t28 = t41 * pkin(3) - t43 * qJ(4) + t78;
t31 = t59 * t47 - t61 * t48;
t12 = t58 * t28 + t60 * t31;
t84 = t58 ^ 2 + t60 ^ 2;
t81 = -0.2e1 * pkin(1) * qJD(2);
t54 = -t61 * pkin(2) - pkin(3);
t5 = t60 * t17 - t58 * t23;
t11 = t60 * t28 - t58 * t31;
t76 = 0.2e1 * t84 * qJD(4);
t75 = t5 * t60 + t6 * t58;
t74 = -t5 * t58 + t6 * t60;
t7 = -t60 * t43 * pkin(7) + t41 * pkin(4) + t11;
t8 = -pkin(7) * t90 + t12;
t73 = t62 * t8 - t64 * t7;
t72 = t62 * t7 + t64 * t8;
t71 = t22 * t43 + t30 * t36;
t70 = t44 * t35 + t41 * t94;
t39 = t92 * t58;
t40 = t92 * t60;
t69 = -t64 * t39 - t62 * t40;
t68 = -t62 * t39 + t64 * t40;
t67 = -qJD(4) * t41 - t35 * t51 + t36 * t54;
t46 = -t60 * pkin(4) + t54;
t25 = t95 * t43;
t24 = t44 * t43;
t21 = pkin(4) * t90 + t30;
t19 = -t44 * qJD(4) - t68 * qJD(5);
t18 = -qJD(4) * t95 - t69 * qJD(5);
t16 = t35 * t95 - t38 * t41;
t13 = pkin(4) * t89 + t22;
t10 = t44 * t36 + t94 * t43;
t9 = t36 * t95 - t43 * t38;
t4 = -pkin(7) * t89 + t6;
t3 = t35 * pkin(4) - pkin(7) * t88 + t5;
t2 = -t72 * qJD(5) + t64 * t3 - t62 * t4;
t1 = t73 * qJD(5) - t62 * t3 - t64 * t4;
t14 = [0, 0, 0, 0.2e1 * t63 * t82, 0.2e1 * (-t63 ^ 2 + t65 ^ 2) * qJD(2), 0, 0, 0, t63 * t81, t65 * t81, -0.2e1 * t23 * t41 - 0.2e1 * t31 * t35 + 0.2e1 * t71, 0.2e1 * t31 * t23 + 0.2e1 * t78 * t55 + 0.2e1 * t91, 0.2e1 * t11 * t35 + 0.2e1 * t5 * t41 + 0.2e1 * t71 * t58, -0.2e1 * t12 * t35 - 0.2e1 * t6 * t41 + 0.2e1 * t71 * t60, -0.2e1 * t75 * t43 + 0.2e1 * (-t11 * t60 - t12 * t58) * t36, 0.2e1 * t11 * t5 + 0.2e1 * t12 * t6 + 0.2e1 * t91, 0.2e1 * t25 * t9, -0.2e1 * t25 * t10 - 0.2e1 * t9 * t24, 0.2e1 * t25 * t35 + 0.2e1 * t9 * t41, -0.2e1 * t10 * t41 - 0.2e1 * t24 * t35, 0.2e1 * t41 * t35, 0.2e1 * t21 * t10 + 0.2e1 * t13 * t24 + 0.2e1 * t2 * t41 - 0.2e1 * t73 * t35, 0.2e1 * t1 * t41 + 0.2e1 * t13 * t25 + 0.2e1 * t21 * t9 - 0.2e1 * t72 * t35; 0, 0, 0, 0, 0, t82, -t83, 0, -pkin(6) * t82, pkin(6) * t83, (-t35 * t59 - t36 * t61) * pkin(2), (-t22 * t61 + t23 * t59) * pkin(2), -t22 * t60 + t67 * t58, t22 * t58 + t67 * t60, t74, t22 * t54 + t74 * t51 + (-t11 * t58 + t12 * t60) * qJD(4), t25 * t94 + t9 * t44, -t44 * t10 - t24 * t94 - t25 * t38 + t9 * t95, t70, t16, 0, t46 * t10 - t13 * t95 + t19 * t41 + t21 * t38 + t69 * t35, t13 * t44 + t18 * t41 + t21 * t94 - t68 * t35 + t46 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t51 * t76, t44 * t93, -0.2e1 * t44 * t38 + 0.2e1 * t94 * t95, 0, 0, 0, 0.2e1 * t46 * t38, t46 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t60 * t35, -t58 * t35, -t84 * t36, t75, 0, 0, 0, 0, 0, t16, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t88, 0, t22, 0, 0, 0, 0, 0, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10, t35, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t38, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
