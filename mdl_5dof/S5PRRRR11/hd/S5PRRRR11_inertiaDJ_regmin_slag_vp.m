% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:45:58
% DurationCPUTime: 0.69s
% Computational Cost: add. (872->100), mult. (2789->177), div. (0->0), fcn. (2545->8), ass. (0->76)
t46 = sin(pkin(5));
t53 = cos(qJ(3));
t80 = qJD(3) * t53;
t65 = t46 * t80;
t85 = t46 * t53;
t92 = qJD(4) * t85 + t65;
t47 = cos(pkin(5));
t91 = -0.2e1 * t47;
t90 = pkin(7) + pkin(8);
t89 = pkin(2) * t47;
t88 = t47 * pkin(4);
t52 = cos(qJ(4));
t87 = t52 * pkin(3);
t49 = sin(qJ(4));
t50 = sin(qJ(3));
t83 = t49 * t50;
t32 = (-t52 * t53 + t83) * t46;
t69 = t46 * t90;
t57 = t50 * t69;
t27 = (pkin(2) * t53 + pkin(3)) * t47 - t57;
t75 = t50 * t89;
t31 = t90 * t85 + t75;
t82 = t52 * t31;
t56 = -t49 * t27 - t82;
t12 = -t32 * pkin(9) - t56;
t51 = cos(qJ(5));
t86 = t12 * t51;
t84 = t49 * t31;
t81 = qJD(3) * t50;
t79 = qJD(4) * t49;
t78 = qJD(4) * t52;
t48 = sin(qJ(5));
t77 = qJD(5) * t48;
t76 = qJD(5) * t51;
t43 = t80 * t89;
t28 = -qJD(3) * t57 + t43;
t29 = (-t53 * t69 - t75) * qJD(3);
t74 = -t27 * t78 - t52 * t28 - t49 * t29;
t73 = pkin(3) * t79;
t72 = pkin(3) * t78;
t71 = pkin(4) * t77;
t70 = pkin(4) * t76;
t45 = t46 ^ 2;
t68 = t45 * t80;
t66 = t46 * t81;
t59 = t46 * t50 * t78 + t92 * t49 + t52 * t66;
t8 = t31 * t79 + t74;
t4 = -t59 * pkin(9) - t8;
t19 = (qJD(3) + qJD(4)) * t46 * t83 - t92 * t52;
t61 = -t49 * t28 + t52 * t29;
t9 = t56 * qJD(4) + t61;
t5 = t19 * pkin(9) + t9;
t64 = -t48 * t4 + t51 * t5;
t33 = (t49 * t53 + t50 * t52) * t46;
t11 = -t33 * pkin(9) + t52 * t27 - t84 + t88;
t63 = -t11 - t88;
t62 = t12 * t77 - t48 * t5;
t44 = pkin(4) + t87;
t60 = qJD(5) * (-pkin(4) - t44);
t58 = pkin(3) * t66;
t36 = (-pkin(3) * t53 - pkin(2)) * t46;
t18 = -t48 * t32 + t51 * t33;
t1 = (-qJD(5) * t11 - t4) * t51 + t62;
t2 = (-t11 * t48 - t86) * qJD(5) + t64;
t55 = (t49 * t77 + (t48 * t49 - t51 * t52) * qJD(4)) * pkin(3);
t54 = (-t49 * t76 + (-t48 * t52 - t49 * t51) * qJD(4)) * pkin(3);
t35 = (-pkin(7) * t85 - t75) * qJD(3);
t34 = pkin(7) * t66 - t43;
t22 = t32 * pkin(4) + t36;
t21 = -t44 * t77 + t54;
t20 = -t44 * t76 + t55;
t17 = t51 * t32 + t48 * t33;
t13 = t59 * pkin(4) + t58;
t7 = t18 * qJD(5) - t48 * t19 + t51 * t59;
t6 = t51 * t19 + t32 * t76 + t33 * t77 + t48 * t59;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0.2e1 * t50 * t68, 0.2e1 * (-t50 ^ 2 + t53 ^ 2) * t45 * qJD(3), 0.2e1 * t47 * t65, t66 * t91, 0, -0.2e1 * t45 * pkin(2) * t81 + 0.2e1 * t35 * t47, -0.2e1 * pkin(2) * t68 + 0.2e1 * t34 * t47, -0.2e1 * t33 * t19, 0.2e1 * t19 * t32 - 0.2e1 * t33 * t59, t19 * t91, t59 * t91, 0, 0.2e1 * t32 * t58 + 0.2e1 * t36 * t59 + 0.2e1 * t9 * t47, -0.2e1 * t36 * t19 + 0.2e1 * t33 * t58 + 0.2e1 * t8 * t47, -0.2e1 * t18 * t6, 0.2e1 * t6 * t17 - 0.2e1 * t18 * t7, t6 * t91, t7 * t91, 0, 0.2e1 * t13 * t17 + 0.2e1 * t2 * t47 + 0.2e1 * t22 * t7, 0.2e1 * t1 * t47 + 0.2e1 * t13 * t18 - 0.2e1 * t22 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, 0, 0, 0, 0, -t59, t19, 0, 0, 0, 0, 0, -t7, t6; 0, 0, 0, 0, 0, 0, t65, -t66, 0, t35, t34, 0, 0, -t19, -t59, 0, (-t82 + (-t47 * pkin(3) - t27) * t49) * qJD(4) + t61, (-t47 * t87 + t84) * qJD(4) + t74, 0, 0, -t6, -t7, 0, t21 * t47 + t2, t20 * t47 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t73, -0.2e1 * t72, 0, 0, 0, 0, 0, 0.2e1 * t21, 0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t19, 0, 0, 0, 0, 0, -t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t59, 0, t9, t8, 0, 0, -t6, -t7, 0, (t63 * t48 - t86) * qJD(5) + t64, (t63 * qJD(5) - t4) * t51 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t72, 0, 0, 0, 0, 0, t48 * t60 + t54, t51 * t60 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t71, -0.2e1 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
