% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPPRR1
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:07
% EndTime: 2019-12-05 14:58:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (136->22), mult. (391->56), div. (0->0), fcn. (320->8), ass. (0->30)
t20 = sin(pkin(9));
t22 = cos(pkin(9));
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t15 = t25 * t20 - t27 * t22;
t11 = t15 * qJD(4);
t44 = t11 * qJD(4);
t16 = t27 * t20 + t25 * t22;
t12 = t16 * qJD(4);
t21 = sin(pkin(8));
t7 = t21 * t12;
t24 = sin(qJ(5));
t26 = cos(qJ(5));
t43 = t24 * t26;
t28 = qJD(5) ^ 2;
t42 = t28 * t24;
t41 = t28 * t26;
t40 = t24 ^ 2 - t26 ^ 2;
t39 = qJD(4) * pkin(4);
t38 = t21 * t44;
t36 = t24 * qJD(5);
t35 = t26 * qJD(5);
t34 = 0.2e1 * qJD(4) * qJD(5);
t33 = t39 * qJD(4);
t32 = pkin(6) * t28;
t31 = -0.2e1 * qJD(5) * t39;
t30 = 0.2e1 * t7 + qJD(5) * cos(pkin(8));
t29 = qJD(4) ^ 2;
t10 = t15 * t21;
t1 = [0, 0, 0, 0, t38, t7 * qJD(4), 0, 0, 0, 0, 0, t26 * t38 + (t10 * t35 + t30 * t24) * qJD(5), -t24 * t38 + (-t10 * t36 + t30 * t26) * qJD(5); 0, 0, 0, 0, -t12 * qJD(4), t44, 0, 0, 0, 0, 0, t11 * t36 - t16 * t41 + (-t12 * t26 + t15 * t36) * qJD(4), t11 * t35 + t16 * t42 + (t12 * t24 + t15 * t35) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t41; 0, 0, 0, 0, 0, 0, t34 * t43, -t40 * t34, t41, -t42, 0, t24 * t31 - t32 * t26, t32 * t24 + t26 * t31; 0, 0, 0, 0, 0, 0, -t29 * t43, t40 * t29, 0, 0, 0, t33 * t24, t33 * t26;];
tauc_reg = t1;
