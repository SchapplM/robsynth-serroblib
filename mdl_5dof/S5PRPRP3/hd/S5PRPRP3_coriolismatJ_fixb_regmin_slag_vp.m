% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x14]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:43
% EndTime: 2019-12-05 15:33:44
% DurationCPUTime: 0.25s
% Computational Cost: add. (342->46), mult. (760->80), div. (0->0), fcn. (743->6), ass. (0->45)
t44 = sin(qJ(4));
t42 = t44 ^ 2;
t45 = cos(qJ(4));
t43 = t45 ^ 2;
t35 = t42 + t43;
t58 = t35 * qJD(2);
t68 = pkin(4) * t44;
t67 = cos(qJ(2));
t66 = sin(qJ(2));
t61 = sin(pkin(8));
t62 = cos(pkin(8));
t30 = t61 * t66 - t62 * t67;
t14 = t30 * t44;
t37 = t61 * t67;
t38 = t62 * t66;
t31 = -t37 - t38;
t3 = (-0.1e1 + t35) * t31 * t30;
t65 = t3 * qJD(1);
t41 = -pkin(2) * t62 - pkin(3);
t32 = -t45 * pkin(4) + t41;
t4 = t32 * t68;
t64 = t4 * qJD(2);
t40 = pkin(2) * t61 + pkin(6);
t63 = qJ(5) + t40;
t60 = qJD(2) * t44;
t59 = qJD(2) * t45;
t36 = t43 - t42;
t57 = t36 * qJD(2);
t56 = t44 * qJD(4);
t55 = t45 * qJD(4);
t54 = pkin(4) * t56;
t53 = pkin(4) * t60;
t52 = t41 * t60;
t51 = t41 * t59;
t50 = t44 * t59;
t49 = t31 * t55;
t48 = t38 / 0.2e1 + t37 / 0.2e1;
t29 = t63 * t45;
t47 = -t29 * t45 - t42 * t63;
t6 = (t43 / 0.2e1 + t42 / 0.2e1) * t31 + t48;
t46 = -qJD(1) * t6 - qJD(2) * t47;
t15 = t30 * t45;
t7 = t48 - t35 * t31 / 0.2e1;
t2 = pkin(4) * t14;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * qJD(2); 0, 0, -t66 * qJD(2), -t67 * qJD(2), (-t30 * t61 + t31 * t62) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, qJD(4) * t14 + t31 * t59, qJD(4) * t15 - t31 * t60, -t30 * t58, t65 + (t30 * t47 - t31 * t32) * qJD(2) + t2 * qJD(4) + t7 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t14 + t49, qJD(2) * t15 - t31 * t56, 0, pkin(4) * t49 + qJD(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t6 - t65; 0, 0, 0, 0, 0, t44 * t55, t36 * qJD(4), 0, 0, 0, t41 * t56, t41 * t55, t35 * qJD(5), qJD(4) * t4 - qJD(5) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t50, t57, t55, -t56, 0, -t40 * t55 + t52, t40 * t56 + t51, -pkin(4) * t55, -pkin(4) * qJD(4) * t29 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t55, 0, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t50, -t57, 0, 0, 0, -t52, -t51, 0, -qJD(5) * t68 - t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t46 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
