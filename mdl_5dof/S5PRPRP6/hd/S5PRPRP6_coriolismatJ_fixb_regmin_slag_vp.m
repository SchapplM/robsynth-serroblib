% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:23
% EndTime: 2019-12-05 15:41:25
% DurationCPUTime: 0.41s
% Computational Cost: add. (164->77), mult. (378->92), div. (0->0), fcn. (283->4), ass. (0->55)
t31 = sin(qJ(4));
t29 = t31 ^ 2;
t33 = cos(qJ(4));
t30 = t33 ^ 2;
t62 = t29 + t30;
t40 = t31 * pkin(4) - t33 * qJ(5);
t66 = t40 * qJD(4) - t31 * qJD(5);
t32 = sin(qJ(2));
t65 = t32 / 0.2e1;
t64 = t33 * pkin(4);
t14 = qJ(3) + t40;
t60 = t31 * qJ(5);
t15 = t60 + t64;
t3 = t14 * t33 + t15 * t31;
t61 = t3 * qJD(2);
t4 = -t14 * t31 + t15 * t33;
t59 = t4 * qJD(2);
t34 = cos(qJ(2));
t5 = (t62 - 0.1e1) * t32 * t34;
t58 = t5 * qJD(1);
t57 = qJD(3) * t33;
t56 = t14 * qJD(2);
t19 = t29 - t30;
t55 = t19 * qJD(2);
t54 = t30 * qJD(2);
t53 = t31 * qJD(4);
t25 = t32 * qJD(2);
t26 = t33 * qJD(2);
t51 = t33 * qJD(4);
t50 = t33 * qJD(5);
t49 = t34 * qJD(2);
t48 = qJ(3) * qJD(4);
t47 = qJD(2) * qJ(3);
t46 = qJD(4) * qJ(5);
t45 = t14 * t26;
t35 = -pkin(2) - pkin(6);
t44 = t35 * t53;
t43 = t35 * t51;
t42 = t31 * t47;
t41 = t33 * t47;
t7 = (0.1e1 / 0.2e1 - t29 / 0.2e1 - t30 / 0.2e1) * t32;
t39 = t7 * qJD(1) + t56;
t37 = t64 / 0.2e1 + t60 / 0.2e1;
t1 = (-t15 / 0.2e1 + t37) * t32;
t36 = t1 * qJD(1) - t15 * t56;
t27 = qJD(3) * t31;
t24 = t31 * qJD(2);
t22 = t31 * t26;
t11 = -t32 * t53 + t33 * t49;
t10 = t31 * t49 + t32 * t51;
t9 = -t31 * t25 + t34 * t51;
t8 = t32 * t26 + t34 * t53;
t6 = (0.1e1 + t62) * t65;
t2 = t15 * t65 + t37 * t32;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * qJD(2); 0, 0, -t25, -t49, t25, t49, (-t32 * pkin(2) + t34 * qJ(3)) * qJD(2) + t32 * qJD(3), 0, 0, 0, 0, 0, t10, t11, t10, -t62 * t25, -t11, t14 * t49 - t58 + t6 * qJD(3) + t2 * qJD(4) + (t62 * qJD(2) * t35 - t50) * t32; 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, t8, 0, -t9, t2 * qJD(2) + t66 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(3) - t1 * qJD(4) + t58; 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), -t31 * t51, t19 * qJD(4), 0, 0, 0, t33 * t48 + t27, -t31 * t48 + t57, t3 * qJD(4) - t31 * t50 + t27, 0, -t4 * qJD(4) + t30 * qJD(5) - t57, (qJD(4) * t15 + qJD(3) - t50) * t14; 0, 0, 0, 0, 0, qJD(2), t47, 0, 0, 0, 0, 0, t24, t26, t24, 0, -t26, t39; 0, 0, 0, 0, 0, 0, 0, -t22, t55, -t53, -t51, 0, t41 - t44, -t42 - t43, -t44 + t61, t66, t43 - t59, -t66 * t35 - t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t53, t54, t44 - t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * qJD(2); 0, 0, 0, 0, 0, -qJD(2), -t47, 0, 0, 0, 0, 0, -t24, -t26, -t24, 0, t26, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t51, -t53, 0, t51, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * qJD(2); 0, 0, 0, 0, 0, 0, 0, t22, -t55, 0, 0, 0, -t41, t42, -t61, 0, t59, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t54, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t12;
