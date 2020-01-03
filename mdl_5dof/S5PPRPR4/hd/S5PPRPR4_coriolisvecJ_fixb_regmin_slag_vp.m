% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:24
% EndTime: 2019-12-31 17:32:25
% DurationCPUTime: 0.25s
% Computational Cost: add. (190->56), mult. (529->93), div. (0->0), fcn. (374->6), ass. (0->44)
t27 = sin(pkin(8));
t28 = cos(pkin(8));
t50 = t27 ^ 2 + t28 ^ 2;
t31 = cos(qJ(5));
t54 = t31 * t28;
t29 = sin(qJ(5));
t55 = t29 * t27;
t13 = -t54 + t55;
t32 = cos(qJ(3));
t43 = t32 * qJD(2);
t38 = qJD(4) - t43;
t30 = sin(qJ(3));
t33 = qJD(3) ^ 2;
t53 = t33 * t30;
t52 = t33 * t32;
t51 = pkin(6) + qJ(4);
t49 = qJD(3) * pkin(3);
t14 = t27 * t31 + t28 * t29;
t9 = qJD(3) * t14;
t48 = qJD(5) * t32;
t10 = t13 * qJD(5);
t47 = t10 * qJD(5);
t11 = t14 * qJD(5);
t46 = t11 * qJD(5);
t45 = t30 * qJD(2);
t44 = t30 * qJD(5) ^ 2;
t42 = qJD(3) * t55;
t41 = qJD(3) * t54;
t24 = -pkin(4) * t28 - pkin(3);
t16 = (qJD(4) + t43) * qJD(3);
t40 = t50 * t16;
t37 = t50 * (qJD(3) * qJ(4) + t45);
t36 = -t37 + t45;
t35 = t14 * t48;
t34 = t13 * t48;
t20 = qJD(5) * t41;
t19 = t38 - t49;
t18 = t51 * t28;
t17 = t51 * t27;
t12 = qJD(3) * t24 + t38;
t7 = -t41 + t42;
t4 = qJD(3) * t11;
t3 = -qJD(5) * t42 + t20;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t47; 0, 0, 0, -t53, -t52, -t28 * t53, t27 * t53, t50 * t52, t30 * t40 + (t19 * t30 - t32 * t36) * qJD(3), 0, 0, 0, 0, 0, -t32 * t4 + t13 * t44 + (t30 * t7 - t35) * qJD(3), -t32 * t3 + t14 * t44 + (t30 * t9 + t34) * qJD(3); 0, 0, 0, 0, 0, 0, 0, qJD(3) * t38 * t50 + t40, t37 * qJD(4) + qJ(4) * t40 + (-t37 * t32 + (-t19 - t49) * t30) * qJD(2), -t10 * t9 + t14 * t3, t10 * t7 - t11 * t9 - t13 * t3 - t14 * t4, -t47, -t46, 0, t12 * t11 + t24 * t4 + ((t17 * t29 - t18 * t31) * qJD(5) - t14 * qJD(4)) * qJD(5) + (t35 + (qJD(3) * t13 - t7) * t30) * qJD(2), -t12 * t10 + t24 * t3 + ((t17 * t31 + t18 * t29) * qJD(5) + t13 * qJD(4)) * qJD(5) - t34 * qJD(2); 0, 0, 0, 0, 0, 0, 0, -t50 * t33, t36 * qJD(3), 0, 0, 0, 0, 0, 0.2e1 * t9 * qJD(5), t20 + (-t7 - t42) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t7, -t7 ^ 2 + t9 ^ 2, t20 + (t7 - t42) * qJD(5), 0, 0, -t12 * t9 - t14 * t16, t12 * t7 + t13 * t16;];
tauc_reg = t1;
