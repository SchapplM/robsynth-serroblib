% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:37
% EndTime: 2019-12-31 17:01:37
% DurationCPUTime: 0.46s
% Computational Cost: add. (279->64), mult. (699->87), div. (0->0), fcn. (551->6), ass. (0->60)
t36 = sin(pkin(7));
t39 = sin(qJ(2));
t29 = t36 * t39 * pkin(1);
t41 = cos(qJ(2));
t63 = t41 * pkin(1);
t32 = pkin(2) + t63;
t37 = cos(pkin(7));
t47 = t37 * t32 - t29;
t16 = -pkin(3) - t47;
t31 = -t37 * pkin(2) - pkin(3);
t71 = t16 / 0.2e1 + t31 / 0.2e1;
t22 = t37 * t63 - t29;
t66 = -t22 / 0.2e1;
t70 = t66 + t71;
t40 = cos(qJ(4));
t51 = qJD(1) + qJD(2);
t69 = t40 * t51;
t38 = sin(qJ(4));
t34 = t38 ^ 2;
t35 = t40 ^ 2;
t27 = t35 - t34;
t68 = t51 * t27;
t62 = t37 * t39;
t61 = pkin(1) * qJD(1);
t60 = pkin(1) * qJD(2);
t45 = pkin(1) * t62 + t36 * t32;
t17 = pkin(6) + t45;
t21 = (t36 * t41 + t62) * pkin(1);
t8 = (t34 + t35) * t22;
t1 = t16 * t21 + t17 * t8;
t59 = t1 * qJD(1);
t2 = -t47 * t21 + t45 * t22;
t58 = t2 * qJD(1);
t57 = t8 * qJD(1);
t56 = qJD(1) * t16;
t55 = qJD(2) * t31;
t54 = t21 * qJD(1);
t18 = t21 * qJD(2);
t53 = t22 * qJD(1);
t52 = t38 * qJD(4);
t33 = t40 * qJD(4);
t50 = t38 * t56;
t49 = t40 * t56;
t48 = t38 * t54;
t46 = pkin(1) * t51;
t44 = t66 - t71;
t3 = t44 * t38;
t43 = t3 * qJD(1) - t38 * t55;
t4 = t44 * t40;
t42 = t4 * qJD(1) - t40 * t55;
t30 = t36 * pkin(2) + pkin(6);
t28 = t38 * t33;
t26 = t27 * qJD(4);
t20 = t38 * t69;
t19 = t22 * qJD(2);
t13 = t38 * t18;
t7 = t8 * qJD(2);
t6 = t70 * t40;
t5 = t70 * t38;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39 * t60, -t41 * t60, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, t2 * qJD(2), t28, t26, 0, -t28, 0, 0, t16 * t52 - t40 * t18, t16 * t33 + t13, t7, t1 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39 * t46, -t41 * t46, 0, 0, 0, 0, 0, 0, 0, 0, -t18 - t54, -t19 - t53, 0, t58 + (-t21 * t37 + t22 * t36) * qJD(2) * pkin(2), t28, t26, 0, -t28, 0, 0, t5 * qJD(4) - t21 * t69, t6 * qJD(4) + t13 + t48, t7 + t57, t59 + (t21 * t31 + t30 * t8) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t68, t33, -t20, -t52, 0, t5 * qJD(2) - t17 * t33 + t50, t6 * qJD(2) + t17 * t52 + t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t61, t41 * t61, 0, 0, 0, 0, 0, 0, 0, 0, t54, t53, 0, -t58, t28, t26, 0, -t28, 0, 0, -t3 * qJD(4) + t40 * t54, -t4 * qJD(4) - t48, -t57, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t26, 0, -t28, 0, 0, t31 * t52, t31 * t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t68, t33, -t20, -t52, 0, -t30 * t33 - t43, t30 * t52 - t42, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t68, 0, t20, 0, 0, t3 * qJD(2) - t50, t4 * qJD(2) - t49, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t68, 0, t20, 0, 0, t43, t42, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
