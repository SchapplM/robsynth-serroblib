% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRPR3
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
% cmat_reg [(5*%NQJ)%x13]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:24
% EndTime: 2019-12-05 15:05:25
% DurationCPUTime: 0.21s
% Computational Cost: add. (226->29), mult. (665->59), div. (0->0), fcn. (762->8), ass. (0->31)
t36 = sin(pkin(9));
t40 = sin(qJ(3));
t54 = cos(pkin(9));
t57 = cos(qJ(3));
t27 = t36 * t40 - t54 * t57;
t55 = qJD(3) * pkin(3);
t39 = sin(qJ(5));
t53 = qJD(3) * t39;
t41 = cos(qJ(5));
t52 = qJD(3) * t41;
t30 = -t39 ^ 2 + t41 ^ 2;
t51 = t30 * qJD(3);
t50 = t39 * qJD(5);
t49 = t40 * qJD(3);
t48 = t41 * qJD(5);
t35 = -t54 * pkin(3) - pkin(4);
t47 = t35 * t53;
t46 = t35 * t52;
t45 = t39 * t52;
t44 = t57 * qJD(3);
t42 = t36 * t57 + t54 * t40;
t38 = cos(pkin(8));
t37 = sin(pkin(8));
t34 = t36 * pkin(3) + pkin(6);
t22 = t42 * t37;
t21 = t27 * t37;
t12 = t27 * t41;
t11 = t27 * t39;
t8 = t22 * t41;
t7 = t22 * t39;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t37 * t44, t37 * t49, (t54 * t21 - t22 * t36) * t55, 0, 0, 0, 0, 0, t7 * qJD(5) + t21 * t52, t8 * qJD(5) - t21 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(3) + (t41 * t21 + t38 * t39) * qJD(5), t8 * qJD(3) + (-t39 * t21 + t38 * t41) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t49, -t44, (-t27 * t36 - t42 * t54) * t55, 0, 0, 0, 0, 0, t11 * qJD(5) - t42 * t52, t12 * qJD(5) + t42 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * qJD(3) - t42 * t48, t12 * qJD(3) + t42 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t39 * t48, t30 * qJD(5), 0, 0, 0, t35 * t50, t35 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t45, t51, t48, -t50, 0, -t34 * t48 + t47, t34 * t50 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t45, -t51, 0, 0, 0, -t47, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
