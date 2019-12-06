% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:44
% EndTime: 2019-12-05 15:24:46
% DurationCPUTime: 0.33s
% Computational Cost: add. (259->56), mult. (627->97), div. (0->0), fcn. (717->8), ass. (0->55)
t40 = cos(pkin(9));
t67 = cos(qJ(5));
t50 = t67 * t40;
t39 = sin(pkin(9));
t41 = sin(qJ(5));
t65 = t41 * t39;
t44 = t50 - t65;
t21 = t44 * qJD(5);
t58 = sin(pkin(8));
t32 = t58 * pkin(2) + qJ(4);
t69 = pkin(6) + t32;
t68 = cos(qJ(2));
t66 = sin(qJ(2));
t64 = t41 * t40;
t37 = t39 ^ 2;
t38 = t40 ^ 2;
t29 = t37 + t38;
t59 = cos(pkin(8));
t23 = t58 * t66 - t59 * t68;
t33 = t58 * t68;
t34 = t59 * t66;
t24 = -t33 - t34;
t1 = (-0.1e1 + t29) * t24 * t23;
t63 = t1 * qJD(1);
t51 = t67 * t39;
t27 = t51 + t64;
t43 = t64 / 0.2e1 + t51 / 0.2e1;
t2 = (-t27 / 0.2e1 + t43) * t23;
t62 = t2 * qJD(1);
t42 = t50 / 0.2e1 - t65 / 0.2e1;
t3 = (-t44 / 0.2e1 + t42) * t23;
t61 = t3 * qJD(1);
t9 = -t27 ^ 2 + t44 ^ 2;
t60 = t9 * qJD(2);
t57 = qJD(2) * t24;
t56 = t44 * qJD(2);
t55 = t27 * qJD(2);
t22 = t27 * qJD(5);
t54 = t29 * qJD(2);
t53 = t44 * t55;
t52 = -t34 / 0.2e1 - t33 / 0.2e1;
t49 = t29 * t23;
t46 = -t59 * pkin(2) - pkin(3);
t28 = -t40 * pkin(4) + t46;
t48 = qJD(2) * t28 + qJD(4);
t47 = (-t38 / 0.2e1 - t37 / 0.2e1) * t24;
t12 = t29 * t32;
t8 = t47 + t52;
t45 = t8 * qJD(1) + t12 * qJD(2);
t20 = t69 * t40;
t19 = t69 * t39;
t7 = t47 - t52;
t5 = (t27 / 0.2e1 + t43) * t23;
t4 = (t44 / 0.2e1 + t42) * t23;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, t1 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, -qJD(2) * t66, -qJD(2) * t68, (-t58 * t23 + t59 * t24) * qJD(2) * pkin(2), t40 * t57, -t39 * t57, -qJD(2) * t49, t63 + (-t24 * t46 - t32 * t49) * qJD(2) + t7 * qJD(4), 0, 0, 0, 0, 0, t5 * qJD(5) + t24 * t56, t4 * qJD(5) - t24 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(2) + t24 * t21, t4 * qJD(2) - t24 * t22; 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(4) - t63, 0, 0, 0, 0, 0, -t2 * qJD(5), -t3 * qJD(5); 0, 0, 0, 0, 0, 0, 0, t29 * qJD(4), t12 * qJD(4), t44 * t22, t9 * qJD(5), 0, 0, 0, t28 * t22, t28 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t54, t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t60, t21, -t22, 0, -t62 + t28 * t55 + (t41 * t19 - t67 * t20) * qJD(5), -t61 + t28 * t56 + (t67 * t19 + t41 * t20) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t21; 0, 0, 0, 0, 0, 0, 0, 0, -t8 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t54, -t45, 0, 0, 0, 0, 0, t22, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(2), t3 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t60, 0, 0, 0, -t27 * t48 + t62, -t44 * t48 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
