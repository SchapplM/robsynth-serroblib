% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% tauc_reg [5x13]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:37
% EndTime: 2019-12-05 14:59:38
% DurationCPUTime: 0.21s
% Computational Cost: add. (132->26), mult. (357->67), div. (0->0), fcn. (276->8), ass. (0->32)
t25 = qJD(4) ^ 2;
t21 = sin(qJ(4));
t24 = qJD(5) ^ 2;
t47 = (t24 + t25) * t21;
t16 = sin(pkin(9));
t23 = cos(qJ(4));
t46 = t16 * t23;
t17 = sin(pkin(8));
t18 = cos(pkin(9));
t45 = t17 * t18;
t20 = sin(qJ(5));
t22 = cos(qJ(5));
t44 = t20 * t22;
t43 = t25 * t21;
t42 = t25 * t23;
t41 = t20 ^ 2 - t22 ^ 2;
t39 = qJD(4) * pkin(4);
t19 = cos(pkin(8));
t10 = -t19 * t21 + t23 * t45;
t38 = t10 * t25;
t37 = (t19 * t23 + t21 * t45) * qJD(4);
t36 = qJD(5) * t10;
t35 = qJD(4) * qJD(5);
t34 = t16 * t42;
t33 = 0.2e1 * t35;
t32 = t39 * qJD(4);
t31 = 0.2e1 * qJD(4) * t16 * t21;
t30 = -0.2e1 * t23 * t35;
t28 = pkin(6) * t24;
t27 = -0.2e1 * qJD(5) * t39;
t26 = -qJD(5) * t16 * t17 + 0.2e1 * t37;
t1 = [0, 0, 0, 0, -t38, t37 * qJD(4), 0, 0, 0, 0, 0, -t22 * t38 + (t26 * t20 - t22 * t36) * qJD(5), t20 * t38 + (t20 * t36 + t26 * t22) * qJD(5); 0, 0, 0, 0, -t34, t16 * t43, 0, 0, 0, 0, 0, -t22 * t34 + (t20 * t31 + (t18 * t20 - t22 * t46) * qJD(5)) * qJD(5), t20 * t34 + (t22 * t31 + (t18 * t22 + t20 * t46) * qJD(5)) * qJD(5); 0, 0, 0, 0, -t43, -t42, 0, 0, 0, 0, 0, t20 * t30 - t22 * t47, t20 * t47 + t22 * t30; 0, 0, 0, 0, 0, 0, t33 * t44, -t41 * t33, t24 * t22, -t24 * t20, 0, t20 * t27 - t28 * t22, t28 * t20 + t22 * t27; 0, 0, 0, 0, 0, 0, -t25 * t44, t41 * t25, 0, 0, 0, t32 * t20, t32 * t22;];
tauc_reg = t1;
