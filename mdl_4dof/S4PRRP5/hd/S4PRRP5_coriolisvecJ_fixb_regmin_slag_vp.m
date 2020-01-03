% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tauc_reg [4x13]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:22
% EndTime: 2019-12-31 16:29:23
% DurationCPUTime: 0.21s
% Computational Cost: add. (161->50), mult. (435->96), div. (0->0), fcn. (227->4), ass. (0->49)
t35 = (qJD(2) * qJD(3));
t53 = -2 * t35;
t19 = sin(qJ(2));
t22 = qJD(3) ^ 2;
t23 = qJD(2) ^ 2;
t52 = (t22 + t23) * t19;
t18 = sin(qJ(3));
t20 = cos(qJ(3));
t38 = t19 * qJD(1);
t12 = qJD(2) * pkin(5) + t38;
t30 = qJ(4) * qJD(2) + t12;
t4 = t30 * t18;
t40 = qJD(3) * pkin(3);
t3 = -t4 + t40;
t5 = t30 * t20;
t27 = t18 * t3 - t20 * t5;
t33 = t18 * t35;
t9 = pkin(3) * t33 + qJD(2) * t38;
t51 = t27 * qJD(2) + t9;
t50 = t3 + t4;
t49 = t22 * t18;
t48 = t22 * t20;
t21 = cos(qJ(2));
t47 = t23 * t21;
t46 = -qJ(4) - pkin(5);
t16 = t18 ^ 2;
t17 = t20 ^ 2;
t45 = t16 - t17;
t44 = t16 + t17;
t42 = qJD(2) * pkin(2);
t34 = -t20 * pkin(3) - pkin(2);
t37 = t21 * qJD(1);
t8 = t34 * qJD(2) + qJD(4) - t37;
t41 = qJD(2) * t8;
t39 = qJD(3) * t18;
t36 = qJ(4) * qJD(3);
t32 = qJD(3) * t46;
t29 = t21 * t53;
t28 = qJD(4) + t37;
t26 = qJD(2) * t42;
t25 = -0.2e1 * qJD(3) * t42;
t1 = -t12 * t39 + (-t18 * t36 + t28 * t20) * qJD(2);
t2 = -t20 * t12 * qJD(3) + (-t28 * t18 - t20 * t36) * qJD(2);
t24 = t1 * t20 - t2 * t18 + (-t18 * t5 - t20 * t3) * qJD(3);
t11 = t46 * t20;
t10 = t46 * t18;
t7 = -t18 * qJD(4) + t20 * t32;
t6 = t20 * qJD(4) + t18 * t32;
t13 = [0, 0, -t23 * t19, -t47, 0, 0, 0, 0, 0, t18 * t29 - t20 * t52, t18 * t52 + t20 * t29, t44 * t47, -t51 * t21 + (t24 + t41) * t19; 0, 0, 0, 0, 0.2e1 * t20 * t33, t45 * t53, t48, -t49, 0, -pkin(5) * t48 + t18 * t25, pkin(5) * t49 + t20 * t25, (-t18 * t7 + t20 * t6 + (-t10 * t20 + t11 * t18) * qJD(3) - t44 * t37) * qJD(2) + t24, -t1 * t11 + t5 * t6 + t2 * t10 + t3 * t7 + t9 * t34 + t8 * pkin(3) * t39 + (-t8 * t19 + t27 * t21) * qJD(1); 0, 0, 0, 0, -t18 * t23 * t20, t45 * t23, 0, 0, 0, t18 * t26, t20 * t26, (-t40 + t50) * t20 * qJD(2), t50 * t5 + (-t18 * t41 + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t23, t51;];
tauc_reg = t13;
