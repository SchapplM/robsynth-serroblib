% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:18
% EndTime: 2019-12-05 17:38:20
% DurationCPUTime: 0.46s
% Computational Cost: add. (295->69), mult. (510->74), div. (0->0), fcn. (516->4), ass. (0->47)
t48 = qJD(4) + qJD(5);
t36 = -pkin(6) + qJ(2);
t69 = -pkin(7) + t36;
t38 = sin(qJ(5));
t39 = sin(qJ(4));
t40 = cos(qJ(4));
t63 = cos(qJ(5));
t24 = t38 * t40 + t63 * t39;
t68 = t48 * t24;
t26 = t69 * t39;
t27 = t69 * t40;
t67 = t48 * (-t63 * t26 - t38 * t27);
t66 = t48 * (t38 * t26 - t63 * t27);
t65 = pkin(4) * t38;
t64 = pkin(4) * t40;
t37 = pkin(1) + qJ(3);
t23 = t38 * t39 - t63 * t40;
t5 = -t23 ^ 2 + t24 ^ 2;
t58 = t5 * qJD(1);
t29 = t39 * pkin(4) + t37;
t6 = -t29 * t23 + t24 * t64;
t57 = t6 * qJD(1);
t7 = -t23 * t64 - t29 * t24;
t56 = t7 * qJD(1);
t55 = qJD(5) * t29;
t54 = t23 * qJD(1);
t53 = t24 * qJD(1);
t28 = t39 ^ 2 - t40 ^ 2;
t52 = t28 * qJD(1);
t51 = t37 * qJD(1);
t31 = t39 * qJD(1);
t50 = t39 * qJD(4);
t32 = t40 * qJD(1);
t49 = t40 * qJD(4);
t47 = t29 * t54;
t46 = t29 * t53;
t45 = t39 * t32;
t44 = t63 * qJD(4);
t43 = t63 * qJD(5);
t41 = -qJD(2) + t51;
t35 = qJ(2) * qJD(1);
t34 = qJ(2) * qJD(2);
t21 = t23 * qJD(2);
t19 = t24 * qJD(2);
t12 = t23 * t53;
t10 = t48 * t23;
t1 = [0, 0, 0, 0, qJD(2), t34, qJD(2), qJD(3), t37 * qJD(3) + t34, -t39 * t49, t28 * qJD(4), 0, 0, 0, qJD(3) * t39 + t37 * t49, qJD(3) * t40 - t37 * t50, t23 * t68, t48 * t5, 0, 0, 0, t24 * qJD(3) + t6 * qJD(4) - t23 * t55, -t23 * qJD(3) + t7 * qJD(4) - t24 * t55; 0, 0, 0, 0, qJD(1), t35, qJD(1), 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, qJD(1), t51, 0, 0, 0, 0, 0, t31, t32, 0, 0, 0, 0, 0, t53, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t52, -t50, -t49, 0, t37 * t32 - t36 * t50, -t37 * t31 - t36 * t49, t12, t58, -t68, t10, 0, t57 + t67, t56 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t58, -t68, t10, 0, -t47 + t67, -t46 + t66; 0, 0, 0, 0, -qJD(1), -t35, -qJD(1), 0, -qJD(3) - t35, 0, 0, 0, 0, 0, -t49, t50, 0, 0, 0, 0, 0, t10, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t31, 0, 0, 0, 0, 0, t54, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t53; 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t41, 0, 0, 0, 0, 0, -t31, -t32, 0, 0, 0, 0, 0, -t53, t54; 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t49, 0, 0, 0, 0, 0, -t68, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t52, 0, 0, 0, -t41 * t40, t41 * t39, -t12, -t58, 0, 0, 0, -t21 - t57, -t19 - t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, 0, 0, 0, 0, -t54, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t65, -pkin(4) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t65, (-t44 - t43) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t58, 0, 0, 0, -t21 + t47, -t19 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t65, pkin(4) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
