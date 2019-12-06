% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:12
% EndTime: 2019-12-05 17:10:14
% DurationCPUTime: 0.40s
% Computational Cost: add. (365->72), mult. (961->117), div. (0->0), fcn. (864->8), ass. (0->74)
t51 = sin(qJ(5));
t52 = sin(qJ(4));
t55 = cos(qJ(5));
t56 = cos(qJ(4));
t30 = t51 * t52 - t55 * t56;
t85 = qJD(4) + qJD(5);
t17 = t85 * t30;
t86 = qJD(2) + qJD(3);
t84 = -pkin(7) - pkin(8);
t57 = cos(qJ(3));
t83 = t57 * pkin(2);
t53 = sin(qJ(3));
t46 = t53 * pkin(2) + pkin(7);
t82 = -pkin(8) - t46;
t81 = cos(qJ(2));
t54 = sin(qJ(2));
t80 = t53 * t54;
t32 = t51 * t56 + t55 * t52;
t18 = t85 * t32;
t72 = t52 * qJD(4);
t67 = pkin(4) * t72;
t74 = pkin(2) * qJD(3);
t69 = t53 * t74;
t36 = t67 + t69;
t48 = -t56 * pkin(4) - pkin(3);
t39 = t48 - t83;
t79 = t39 * t18 + t36 * t30;
t78 = -t39 * t17 + t36 * t32;
t77 = t48 * t18 + t30 * t67;
t76 = -t48 * t17 + t32 * t67;
t47 = -pkin(3) - t83;
t49 = t56 * qJD(4);
t75 = t47 * t49 + t52 * t69;
t73 = pkin(4) * qJD(5);
t71 = pkin(3) * t72;
t70 = pkin(3) * t49;
t68 = t57 * t74;
t66 = t51 * t73;
t65 = t55 * t73;
t64 = t57 * t81;
t63 = qJD(4) * t84;
t62 = qJD(4) * t82;
t61 = t81 * qJD(2);
t60 = t52 * t68;
t59 = t56 * t68;
t58 = t47 * t72 - t56 * t69;
t33 = t53 * t81 + t57 * t54;
t50 = t56 * pkin(8);
t43 = 0.2e1 * t52 * t49;
t41 = t56 * pkin(7) + t50;
t40 = t84 * t52;
t35 = t56 * t63;
t34 = t52 * t63;
t31 = -t64 + t80;
t29 = 0.2e1 * (-t52 ^ 2 + t56 ^ 2) * qJD(4);
t28 = t56 * t46 + t50;
t27 = t82 * t52;
t24 = t56 * t62 - t60;
t23 = t52 * t62 + t59;
t20 = t86 * t33;
t19 = -qJD(3) * t64 - t57 * t61 + t86 * t80;
t12 = -0.2e1 * t32 * t17;
t11 = -t20 * t56 + t31 * t72;
t10 = t20 * t52 + t31 * t49;
t9 = -t51 * t34 + t55 * t35 + (-t40 * t51 - t41 * t55) * qJD(5);
t8 = -t55 * t34 - t51 * t35 + (-t40 * t55 + t41 * t51) * qJD(5);
t7 = -t31 * t17 + t20 * t32;
t6 = t31 * t18 + t20 * t30;
t5 = 0.2e1 * t17 * t30 - 0.2e1 * t32 * t18;
t4 = -t51 * t23 + t55 * t24 + (-t27 * t51 - t28 * t55) * qJD(5);
t3 = -t55 * t23 - t51 * t24 + (-t27 * t55 + t28 * t51) * qJD(5);
t2 = t33 * t17 + t32 * t19;
t1 = t18 * t33 - t30 * t19;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t54 * qJD(2), -t61, 0, -t20, t19, 0, 0, 0, 0, 0, t11, t10, 0, 0, 0, 0, 0, t6, t7; 0, 0, 0, 0, 0, -0.2e1 * t69, -0.2e1 * t68, t43, t29, 0, 0, 0, 0.2e1 * t58, 0.2e1 * t75, t12, t5, 0, 0, 0, 0.2e1 * t79, 0.2e1 * t78; 0, 0, 0, 0, 0, -t20, t19, 0, 0, 0, 0, 0, t11, t10, 0, 0, 0, 0, 0, t6, t7; 0, 0, 0, 0, 0, -t69, -t68, t43, t29, 0, 0, 0, t58 - t71, -t70 + t75, t12, t5, 0, 0, 0, t77 + t79, t76 + t78; 0, 0, 0, 0, 0, 0, 0, t43, t29, 0, 0, 0, -0.2e1 * t71, -0.2e1 * t70, t12, t5, 0, 0, 0, 0.2e1 * t77, 0.2e1 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t19 - t33 * t49, t56 * t19 + t33 * t72, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t72, 0, -t46 * t49 - t60, t46 * t72 - t59, 0, 0, -t17, -t18, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t72, 0, -pkin(7) * t49, pkin(7) * t72, 0, 0, -t17, -t18, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t66, -0.2e1 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
