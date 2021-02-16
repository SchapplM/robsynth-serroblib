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
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 19:47:32
% EndTime: 2021-01-15 19:47:36
% DurationCPUTime: 0.75s
% Computational Cost: add. (993->109), mult. (2406->231), div. (0->0), fcn. (2222->8), ass. (0->75)
t59 = sin(pkin(9));
t61 = cos(pkin(9));
t63 = sin(qJ(5));
t65 = cos(qJ(5));
t95 = -t63 * t59 + t65 * t61;
t94 = t95 * qJD(5);
t44 = t59 * t65 + t61 * t63;
t38 = t44 * qJD(5);
t93 = 0.2e1 * t94;
t60 = sin(pkin(8));
t51 = pkin(2) * t60 + qJ(4);
t92 = pkin(7) + t51;
t64 = sin(qJ(2));
t66 = cos(qJ(2));
t85 = -qJ(3) - pkin(6);
t78 = qJD(2) * t85;
t34 = t66 * qJD(3) + t64 * t78;
t62 = cos(pkin(8));
t67 = -t64 * qJD(3) + t66 * t78;
t22 = t34 * t60 - t62 * t67;
t47 = t85 * t64;
t48 = t85 * t66;
t30 = -t47 * t62 - t48 * t60;
t91 = t30 * t22;
t43 = t60 * t66 + t62 * t64;
t90 = t43 * t59;
t82 = t66 * qJD(2);
t83 = t64 * qJD(2);
t36 = -t60 * t83 + t62 * t82;
t89 = t59 * t36;
t88 = t61 * t36;
t35 = t43 * qJD(2);
t56 = pkin(2) * t83;
t17 = pkin(3) * t35 - qJ(4) * t36 - qJD(4) * t43 + t56;
t23 = t62 * t34 + t60 * t67;
t6 = t17 * t59 + t23 * t61;
t41 = t60 * t64 - t62 * t66;
t55 = -pkin(2) * t66 - pkin(1);
t28 = pkin(3) * t41 - qJ(4) * t43 + t55;
t31 = t47 * t60 - t48 * t62;
t12 = t28 * t59 + t31 * t61;
t84 = t59 ^ 2 + t61 ^ 2;
t81 = -0.2e1 * pkin(1) * qJD(2);
t54 = -pkin(2) * t62 - pkin(3);
t5 = t17 * t61 - t23 * t59;
t11 = t28 * t61 - t31 * t59;
t77 = 0.2e1 * t84 * qJD(4);
t76 = t5 * t61 + t59 * t6;
t75 = -t5 * t59 + t6 * t61;
t7 = -pkin(7) * t43 * t61 + pkin(4) * t41 + t11;
t8 = -pkin(7) * t90 + t12;
t74 = t63 * t8 - t65 * t7;
t73 = t63 * t7 + t65 * t8;
t72 = t22 * t43 + t30 * t36;
t71 = t35 * t44 + t41 * t94;
t39 = t92 * t59;
t40 = t92 * t61;
t70 = -t39 * t65 - t40 * t63;
t69 = -t39 * t63 + t40 * t65;
t68 = -qJD(4) * t41 - t35 * t51 + t36 * t54;
t46 = -pkin(4) * t61 + t54;
t25 = t95 * t43;
t24 = t44 * t43;
t21 = pkin(4) * t90 + t30;
t19 = -qJD(4) * t44 - qJD(5) * t69;
t18 = -qJD(4) * t95 - qJD(5) * t70;
t16 = t35 * t95 - t38 * t41;
t13 = pkin(4) * t89 + t22;
t10 = t44 * t36 + t43 * t94;
t9 = t36 * t95 - t38 * t43;
t4 = -pkin(7) * t89 + t6;
t3 = pkin(4) * t35 - pkin(7) * t88 + t5;
t2 = -qJD(5) * t73 + t65 * t3 - t63 * t4;
t1 = qJD(5) * t74 - t63 * t3 - t65 * t4;
t14 = [0, 0, 0, 0.2e1 * t64 * t82, 0.2e1 * (-t64 ^ 2 + t66 ^ 2) * qJD(2), 0, 0, 0, t64 * t81, t66 * t81, 0.2e1 * t35 * t55 + 0.2e1 * t41 * t56, 0.2e1 * t36 * t55 + 0.2e1 * t43 * t56, -0.2e1 * t23 * t41 - 0.2e1 * t31 * t35 + 0.2e1 * t72, 0.2e1 * t23 * t31 + 0.2e1 * t55 * t56 + 0.2e1 * t91, 0.2e1 * t11 * t35 + 0.2e1 * t5 * t41 + 0.2e1 * t59 * t72, -0.2e1 * t12 * t35 - 0.2e1 * t6 * t41 + 0.2e1 * t61 * t72, -0.2e1 * t76 * t43 + 0.2e1 * (-t11 * t61 - t12 * t59) * t36, 0.2e1 * t11 * t5 + 0.2e1 * t12 * t6 + 0.2e1 * t91, 0.2e1 * t25 * t9, -0.2e1 * t10 * t25 - 0.2e1 * t24 * t9, 0.2e1 * t25 * t35 + 0.2e1 * t41 * t9, -0.2e1 * t10 * t41 - 0.2e1 * t24 * t35, 0.2e1 * t41 * t35, 0.2e1 * t10 * t21 + 0.2e1 * t13 * t24 + 0.2e1 * t2 * t41 - 0.2e1 * t35 * t74, 0.2e1 * t1 * t41 + 0.2e1 * t13 * t25 + 0.2e1 * t21 * t9 - 0.2e1 * t35 * t73; 0, 0, 0, 0, 0, t82, -t83, 0, -pkin(6) * t82, pkin(6) * t83, -t22, -t23, (-t35 * t60 - t36 * t62) * pkin(2), (-t22 * t62 + t23 * t60) * pkin(2), -t22 * t61 + t59 * t68, t22 * t59 + t61 * t68, t75, t22 * t54 + t75 * t51 + (-t11 * t59 + t12 * t61) * qJD(4), t25 * t94 + t44 * t9, -t10 * t44 - t24 * t94 - t25 * t38 + t9 * t95, t71, t16, 0, t10 * t46 - t13 * t95 + t19 * t41 + t21 * t38 + t35 * t70, t13 * t44 + t18 * t41 + t21 * t94 - t35 * t69 + t46 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t51 * t77, t44 * t93, -0.2e1 * t38 * t44 + 0.2e1 * t94 * t95, 0, 0, 0, 0.2e1 * t46 * t38, t46 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t36, 0, t56, t61 * t35, -t59 * t35, -t84 * t36, t76, 0, 0, 0, 0, 0, t16, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t88, 0, t22, 0, 0, 0, 0, 0, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10, t35, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t38, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
