% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:26
% EndTime: 2019-12-31 17:32:26
% DurationCPUTime: 0.24s
% Computational Cost: add. (136->51), mult. (390->85), div. (0->0), fcn. (385->6), ass. (0->45)
t31 = sin(pkin(8));
t29 = t31 ^ 2;
t32 = cos(pkin(8));
t30 = t32 ^ 2;
t23 = t29 + t30;
t55 = cos(qJ(5));
t33 = sin(qJ(5));
t54 = t33 * t31;
t53 = t33 * t32;
t52 = pkin(6) + qJ(4);
t40 = t55 * t32;
t13 = -t40 + t54;
t41 = t55 * t31;
t15 = t41 + t53;
t1 = t13 ^ 2 - t15 ^ 2;
t51 = t1 * qJD(3);
t35 = cos(qJ(3));
t37 = -t53 / 0.2e1 - t41 / 0.2e1;
t2 = (t15 / 0.2e1 + t37) * t35;
t50 = t2 * qJD(2);
t36 = -t40 / 0.2e1 + t54 / 0.2e1;
t3 = (-t13 / 0.2e1 + t36) * t35;
t49 = t3 * qJD(2);
t34 = sin(qJ(3));
t6 = (-0.1e1 + t23) * t35 * t34;
t48 = t6 * qJD(2);
t47 = qJD(3) * t34;
t46 = qJD(3) * t35;
t45 = t13 * qJD(3);
t11 = t13 * qJD(5);
t44 = t15 * qJD(3);
t12 = t15 * qJD(5);
t43 = t23 * qJD(3);
t42 = t13 * t44;
t16 = t23 * qJ(4);
t27 = -t32 * pkin(4) - pkin(3);
t39 = qJD(3) * t27 + qJD(4);
t7 = (0.1e1 / 0.2e1 - t30 / 0.2e1 - t29 / 0.2e1) * t34;
t38 = t7 * qJD(2) - t16 * qJD(3);
t22 = t52 * t32;
t21 = t52 * t31;
t8 = (0.1e1 + t23) * t34 / 0.2e1;
t5 = (-t15 / 0.2e1 + t37) * t35;
t4 = (t13 / 0.2e1 + t36) * t35;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t6 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t47, -t46, -t32 * t47, t31 * t47, t23 * t46, t48 + (-t34 * pkin(3) + t16 * t35) * qJD(3) + t8 * qJD(4), 0, 0, 0, 0, 0, t5 * qJD(5) + t34 * t45, t4 * qJD(5) + t34 * t44; 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(3) + t11 * t34, t4 * qJD(3) + t12 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t7 * qJD(4) - t48, 0, 0, 0, 0, 0, -t2 * qJD(5), -t3 * qJD(5); 0, 0, 0, 0, 0, 0, 0, t23 * qJD(4), t16 * qJD(4), -t13 * t12, t1 * qJD(5), 0, 0, 0, t27 * t12, -t27 * t11; 0, 0, 0, 0, 0, 0, 0, t43, -t38, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t51, -t11, -t12, 0, -t50 + t27 * t44 + (t33 * t21 - t22 * t55) * qJD(5), -t49 - t27 * t45 + (t21 * t55 + t33 * t22) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t43, t38, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(3), t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t51, 0, 0, 0, -t15 * t39 + t50, t13 * t39 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
