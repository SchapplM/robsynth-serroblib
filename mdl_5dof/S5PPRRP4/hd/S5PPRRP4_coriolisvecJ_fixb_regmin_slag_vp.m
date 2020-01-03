% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x14]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:41
% EndTime: 2019-12-31 17:34:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (220->58), mult. (542->104), div. (0->0), fcn. (294->4), ass. (0->52)
t38 = (qJD(3) * qJD(4));
t56 = -2 * t38;
t22 = sin(qJ(3));
t25 = qJD(4) ^ 2;
t26 = qJD(3) ^ 2;
t55 = (t25 + t26) * t22;
t21 = sin(qJ(4));
t23 = cos(qJ(4));
t42 = t22 * qJD(2);
t12 = qJD(3) * pkin(6) + t42;
t33 = qJ(5) * qJD(3) + t12;
t41 = t23 * qJD(1);
t4 = -t33 * t21 - t41;
t44 = qJD(4) * pkin(4);
t3 = t4 + t44;
t17 = t21 * qJD(1);
t5 = t33 * t23 - t17;
t30 = t21 * t3 - t23 * t5;
t36 = t21 * t38;
t9 = pkin(4) * t36 + qJD(3) * t42;
t54 = t30 * qJD(3) + t9;
t53 = t3 - t4;
t52 = t25 * t21;
t18 = t25 * t23;
t24 = cos(qJ(3));
t51 = t26 * t24;
t50 = -qJ(5) - pkin(6);
t19 = t21 ^ 2;
t20 = t23 ^ 2;
t49 = t19 - t20;
t48 = t19 + t20;
t46 = qJD(3) * pkin(3);
t37 = -t23 * pkin(4) - pkin(3);
t40 = t24 * qJD(2);
t8 = t37 * qJD(3) + qJD(5) - t40;
t45 = qJD(3) * t8;
t43 = t21 * qJD(4);
t39 = qJ(5) * qJD(4);
t35 = qJD(4) * t50;
t32 = t24 * t56;
t31 = qJD(5) + t40;
t29 = t46 * qJD(3);
t28 = -0.2e1 * qJD(4) * t46;
t1 = (-t21 * t12 - t41) * qJD(4) + (-t21 * t39 + t31 * t23) * qJD(3);
t16 = qJD(1) * t43;
t2 = -qJD(4) * t23 * t12 + t16 + (-t31 * t21 - t23 * t39) * qJD(3);
t27 = t1 * t23 - t2 * t21 + (-t21 * t5 - t23 * t3) * qJD(4);
t11 = t50 * t23;
t10 = t50 * t21;
t7 = -t21 * qJD(5) + t23 * t35;
t6 = t23 * qJD(5) + t21 * t35;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t18, 0, t30 * qJD(4) - t1 * t21 - t2 * t23; 0, 0, 0, -t26 * t22, -t51, 0, 0, 0, 0, 0, t21 * t32 - t23 * t55, t21 * t55 + t23 * t32, t48 * t51, -t54 * t24 + (t27 + t45) * t22; 0, 0, 0, 0, 0, 0.2e1 * t23 * t36, t49 * t56, t18, -t52, 0, -pkin(6) * t18 + t21 * t28, pkin(6) * t52 + t23 * t28, (-t21 * t7 + t23 * t6 + (-t10 * t23 + t11 * t21) * qJD(4) - t48 * t40) * qJD(3) + t27, -t1 * t11 + t5 * t6 + t2 * t10 + t3 * t7 + t9 * t37 + t8 * pkin(4) * t43 + (-t8 * t22 + t30 * t24) * qJD(2); 0, 0, 0, 0, 0, -t21 * t26 * t23, t49 * t26, 0, 0, 0, -t17 * qJD(4) + t21 * t29 + t16, t23 * t29, (-t44 + t53) * t23 * qJD(3), t53 * t5 + (-t21 * t45 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t26, t54;];
tauc_reg = t13;
