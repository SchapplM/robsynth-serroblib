% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:50
% EndTime: 2019-03-08 18:31:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (284->41), mult. (556->55), div. (0->0), fcn. (452->6), ass. (0->43)
t54 = -pkin(3) / 0.2e1;
t53 = pkin(1) * sin(pkin(7));
t52 = cos(qJ(3));
t32 = sin(qJ(3));
t34 = cos(pkin(7)) * pkin(1) + pkin(2);
t25 = t32 * t53 - t34 * t52;
t24 = pkin(3) - t25;
t31 = sin(qJ(4));
t51 = t31 * t24;
t50 = t31 * t25;
t26 = t32 * t34 + t52 * t53;
t49 = t31 * t26;
t33 = cos(qJ(4));
t48 = t33 * t24;
t47 = t33 * t25;
t21 = t33 * t26;
t46 = pkin(3) * qJD(3);
t45 = pkin(3) * qJD(4);
t35 = t25 / 0.2e1 + pkin(3) / 0.2e1 + t24 / 0.2e1;
t1 = t35 * t31;
t44 = t1 * qJD(1);
t3 = t35 * t33;
t43 = t3 * qJD(1);
t9 = -t48 + t49;
t42 = t9 * qJD(1);
t10 = t21 + t51;
t41 = t10 * qJD(1);
t11 = t21 - t50;
t40 = t11 * qJD(1);
t12 = -t47 - t49;
t39 = t12 * qJD(1);
t38 = t25 * qJD(1);
t37 = t26 * qJD(1);
t36 = pkin(3) * (-qJD(3) - qJD(4));
t23 = t26 * qJD(3);
t22 = t25 * qJD(3);
t8 = t12 * qJD(3);
t7 = t11 * qJD(3);
t6 = t10 * qJD(4);
t5 = t9 * qJD(4);
t4 = t33 * t54 + t49 - t48 / 0.2e1 + t47 / 0.2e1;
t2 = t31 * t54 - t21 - t51 / 0.2e1 + t50 / 0.2e1;
t13 = [0, 0, 0, 0, 0, -t23, t22, 0, -t7 - t6, -t8 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t23 - t37, t22 + t38, 0, qJD(4) * t2 - t40 - t7, qJD(4) * t4 - t39 - t8; 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t2 - t41 - t6, qJD(3) * t4 + t42 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t37, -t38, 0, -qJD(4) * t1 + t40, -qJD(4) * t3 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t31 * t45, -t33 * t45; 0, 0, 0, 0, 0, 0, 0, 0, t31 * t36 - t44, t33 * t36 - t43; 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t1 + t41, qJD(3) * t3 - t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t31 * t46 + t44, t33 * t46 + t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t13;
