% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [5x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:22
% EndTime: 2019-12-05 15:03:22
% DurationCPUTime: 0.15s
% Computational Cost: add. (134->44), mult. (376->69), div. (0->0), fcn. (262->6), ass. (0->44)
t39 = qJD(3) * qJ(4);
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t23 = sin(qJ(3));
t25 = cos(qJ(3));
t11 = t25 * t20 + t23 * t21;
t7 = t11 * qJD(1);
t4 = t7 + t39;
t38 = qJD(1) * qJD(3);
t34 = t25 * t38;
t35 = t23 * t38;
t5 = t20 * t34 + t21 * t35;
t33 = t4 * qJD(3) - t5;
t50 = t25 * t21;
t51 = t23 * t20;
t10 = -t50 + t51;
t14 = t21 * t34;
t36 = qJD(1) * t51;
t6 = -qJD(1) * t50 + t36;
t53 = (-t6 + t36) * qJD(3) - t14;
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t52 = t22 * t24;
t27 = qJD(5) ^ 2;
t49 = t27 * t22;
t48 = t27 * t24;
t47 = t22 ^ 2 - t24 ^ 2;
t28 = qJD(3) ^ 2;
t46 = -t27 - t28;
t8 = t10 * qJD(3);
t44 = t8 * qJD(3);
t9 = t11 * qJD(3);
t43 = t9 * qJD(3);
t42 = qJD(4) + t6;
t41 = t22 * qJD(5);
t40 = t24 * qJD(5);
t37 = qJD(3) * qJD(5);
t31 = t4 - t7 + t39;
t30 = t7 * qJD(3) - t5;
t17 = qJD(3) * qJD(4);
t2 = t20 * t35 - t14 - t17;
t29 = t42 * qJD(3) - (-pkin(3) - pkin(6)) * t27 - t2;
t3 = -qJD(3) * pkin(3) + t42;
t1 = [0, 0, 0, -t43, t44, t43, -t44, t5 * t10 - t2 * t11 + t3 * t9 - t4 * t8, 0, 0, 0, 0, 0, t9 * t40 - t10 * t49 + (t11 * t40 - t22 * t8) * qJD(3), -t9 * t41 - t10 * t48 + (-t11 * t41 - t24 * t8) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t49; 0, 0, 0, t30, t53, -t30, 0.2e1 * t17 - t53, -t5 * pkin(3) - t2 * qJ(4) - t3 * t7 + t42 * t4, -0.2e1 * t37 * t52, 0.2e1 * t47 * t37, -t49, -t48, 0, t29 * t22 + t31 * t40, t29 * t24 - t31 * t41; 0, 0, 0, 0, 0, 0, -t28, -t33, 0, 0, 0, 0, 0, t46 * t22, t46 * t24; 0, 0, 0, 0, 0, 0, 0, 0, t28 * t52, -t47 * t28, 0, 0, 0, -t33 * t24, t33 * t22;];
tauc_reg = t1;
