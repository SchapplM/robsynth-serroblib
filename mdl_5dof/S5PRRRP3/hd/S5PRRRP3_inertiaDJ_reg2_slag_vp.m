% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:13
% EndTime: 2019-12-05 16:44:16
% DurationCPUTime: 0.61s
% Computational Cost: add. (628->86), mult. (1557->148), div. (0->0), fcn. (1332->4), ass. (0->55)
t43 = sin(qJ(4));
t44 = sin(qJ(3));
t65 = t43 * t44;
t72 = qJD(3) + qJD(4);
t73 = t72 * t65;
t71 = -pkin(7) - pkin(6);
t70 = pkin(3) * t43;
t45 = cos(qJ(3));
t68 = cos(qJ(4));
t31 = t43 * t45 + t68 * t44;
t19 = t72 * t31;
t69 = t19 * pkin(4);
t50 = t68 * t45;
t30 = -t50 + t65;
t67 = t30 * t19;
t49 = t68 * qJD(4);
t18 = -qJD(3) * t50 - t45 * t49 + t73;
t66 = t31 * t18;
t48 = pkin(3) * t49;
t63 = -t19 * t70 - t30 * t48;
t62 = -t18 * t70 + t31 * t48;
t35 = t71 * t45;
t55 = t43 * t71;
t21 = -t68 * t35 + t44 * t55;
t61 = qJD(4) * t43;
t60 = t44 * qJD(3);
t59 = t45 * qJD(3);
t58 = -0.2e1 * pkin(2) * qJD(3);
t57 = pkin(3) * t60;
t56 = pkin(3) * t61;
t54 = t68 * pkin(3);
t53 = t30 * t61;
t52 = t31 * t61;
t51 = t44 * t59;
t42 = -t45 * pkin(3) - pkin(2);
t47 = t71 * t68;
t33 = t44 * t47;
t20 = t43 * t35 + t33;
t46 = qJD(3) * t47;
t6 = -qJD(4) * t33 - t35 * t61 - t44 * t46 - t55 * t59;
t7 = t35 * t49 + t45 * t46 - t71 * t73;
t41 = t54 + pkin(4);
t39 = -0.2e1 * t48;
t38 = -0.2e1 * t56;
t22 = t30 * pkin(4) + t42;
t14 = t57 + t69;
t13 = -t30 * qJ(5) + t21;
t12 = -t31 * qJ(5) + t20;
t9 = -0.2e1 * t66;
t8 = 0.2e1 * t67;
t5 = -0.2e1 * t66 + 0.2e1 * t67;
t4 = 0.2e1 * t30 * t18 - 0.2e1 * t31 * t19;
t3 = t18 * qJ(5) - t31 * qJD(5) + t7;
t2 = t19 * qJ(5) + t30 * qJD(5) + t6;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18 * t21 - t19 * t20 - t30 * t7 - t31 * t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t12 - t18 * t13 - t31 * t2 - t30 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t51, 0.2e1 * (-t44 ^ 2 + t45 ^ 2) * qJD(3), 0, -0.2e1 * t51, 0, 0, t44 * t58, t45 * t58, 0, 0, t9, t4, 0, t8, 0, 0, 0.2e1 * t42 * t19 + 0.2e1 * t30 * t57, -0.2e1 * t42 * t18 + 0.2e1 * t31 * t57, 0.2e1 * t20 * t18 - 0.2e1 * t21 * t19 + 0.2e1 * t6 * t30 - 0.2e1 * t7 * t31, 0.2e1 * t20 * t7 - 0.2e1 * t21 * t6 + 0.2e1 * t42 * t57, t9, t4, 0, t8, 0, 0, 0.2e1 * t14 * t30 + 0.2e1 * t22 * t19, 0.2e1 * t14 * t31 - 0.2e1 * t22 * t18, 0.2e1 * t12 * t18 - 0.2e1 * t13 * t19 + 0.2e1 * t2 * t30 - 0.2e1 * t3 * t31, 0.2e1 * t12 * t3 - 0.2e1 * t13 * t2 + 0.2e1 * t22 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t59, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t18, 0, (-t68 * t19 + t53) * pkin(3) + t62, 0, 0, 0, 0, 0, 0, -t19, t18, 0, pkin(3) * t53 - t19 * t41 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, -t60, 0, -pkin(6) * t59, pkin(6) * t60, 0, 0, 0, 0, -t18, 0, -t19, 0, t7, t6, (t68 * t18 + t52) * pkin(3) + t63, (t68 * t7 - t43 * t6 + (-t20 * t43 + t68 * t21) * qJD(4)) * pkin(3), 0, 0, -t18, 0, -t19, 0, t3, t2, pkin(3) * t52 + t41 * t18 + t63, t3 * t41 + (-t2 * t43 + (-t12 * t43 + t68 * t13) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t39, 0, 0, 0, 0, 0, 0, 0, 0, t38, t39, 0, 0.2e1 * (t54 - t41) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t18, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t18, 0, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, -t19, 0, t7, t6, 0, 0, 0, 0, -t18, 0, -t19, 0, t3, t2, t18 * pkin(4), t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t48, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t48, 0, -pkin(4) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
