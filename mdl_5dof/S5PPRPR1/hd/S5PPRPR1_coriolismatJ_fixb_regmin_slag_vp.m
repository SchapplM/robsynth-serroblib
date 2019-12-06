% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:23
% EndTime: 2019-12-05 15:01:25
% DurationCPUTime: 0.27s
% Computational Cost: add. (221->55), mult. (570->90), div. (0->0), fcn. (646->8), ass. (0->53)
t65 = cos(qJ(3));
t64 = cos(qJ(5));
t35 = sin(pkin(9));
t38 = sin(qJ(5));
t63 = t38 * t35;
t37 = cos(pkin(9));
t62 = t38 * t37;
t61 = pkin(6) + qJ(4);
t33 = t35 ^ 2;
t34 = t37 ^ 2;
t27 = t33 + t34;
t36 = sin(pkin(8));
t39 = sin(qJ(3));
t56 = cos(pkin(8));
t21 = t39 * t36 - t65 * t56;
t46 = t39 * t56;
t50 = t65 * t36;
t23 = t50 + t46;
t1 = (0.1e1 - t27) * t23 * t21;
t60 = t1 * qJD(1);
t49 = t64 * t35;
t22 = t49 + t62;
t42 = t62 / 0.2e1 + t49 / 0.2e1;
t2 = (-t22 / 0.2e1 + t42) * t21;
t59 = t2 * qJD(1);
t48 = t64 * t37;
t19 = -t48 + t63;
t41 = t48 / 0.2e1 - t63 / 0.2e1;
t3 = (t19 / 0.2e1 + t41) * t21;
t58 = t3 * qJD(1);
t8 = t19 ^ 2 - t22 ^ 2;
t57 = t8 * qJD(3);
t55 = t19 * qJD(3);
t17 = t19 * qJD(5);
t54 = t22 * qJD(3);
t18 = t22 * qJD(5);
t53 = t23 * qJD(3);
t52 = t27 * qJD(3);
t51 = t19 * t54;
t47 = t27 * t21;
t32 = -t37 * pkin(4) - pkin(3);
t45 = qJD(3) * t32 + qJD(4);
t44 = (t34 / 0.2e1 + t33 / 0.2e1) * t23;
t24 = t27 * qJ(4);
t40 = -t50 / 0.2e1 - t46 / 0.2e1;
t7 = t44 + t40;
t43 = t7 * qJD(1) + t24 * qJD(3);
t26 = t61 * t37;
t25 = t61 * t35;
t6 = t44 - t40;
t5 = (t22 / 0.2e1 + t42) * t21;
t4 = (-t19 / 0.2e1 + t41) * t21;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, t1 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t53, t21 * qJD(3), -t37 * t53, t35 * t53, -qJD(3) * t47, t60 + (-t23 * pkin(3) - qJ(4) * t47) * qJD(3) + t6 * qJD(4), 0, 0, 0, 0, 0, t5 * qJD(5) + t19 * t53, t4 * qJD(5) + t22 * t53; 0, 0, 0, 0, 0, 0, 0, 0, t6 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(3) + t23 * t17, t4 * qJD(3) + t23 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(4) - t60, 0, 0, 0, 0, 0, -t2 * qJD(5), -t3 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t27 * qJD(4), t24 * qJD(4), -t19 * t18, t8 * qJD(5), 0, 0, 0, t32 * t18, -t32 * t17; 0, 0, 0, 0, 0, 0, 0, t52, t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t57, -t17, -t18, 0, -t59 + t32 * t54 + (t38 * t25 - t64 * t26) * qJD(5), -t58 - t32 * t55 + (t64 * t25 + t38 * t26) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, -t7 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t52, -t43, 0, 0, 0, 0, 0, t18, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(3), t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t57, 0, 0, 0, -t45 * t22 + t59, t45 * t19 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
