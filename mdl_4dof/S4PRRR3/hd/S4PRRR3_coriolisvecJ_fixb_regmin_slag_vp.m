% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x14]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:39
% EndTime: 2019-12-31 16:31:39
% DurationCPUTime: 0.14s
% Computational Cost: add. (89->29), mult. (209->54), div. (0->0), fcn. (90->4), ass. (0->34)
t32 = pkin(2) * qJD(3);
t24 = qJD(2) * t32;
t17 = cos(qJ(4));
t30 = t17 * qJD(4);
t15 = sin(qJ(4));
t16 = sin(qJ(3));
t37 = t15 * t16;
t12 = qJD(2) + qJD(3);
t18 = cos(qJ(3));
t33 = pkin(2) * qJD(2);
t6 = -t12 * pkin(3) - t18 * t33;
t38 = t24 * t37 + t6 * t30;
t36 = t16 * t17;
t19 = qJD(4) ^ 2;
t35 = t19 * t15;
t10 = t19 * t17;
t34 = t15 ^ 2 - t17 ^ 2;
t31 = t15 * qJD(4);
t29 = -qJD(2) - t12;
t28 = -qJD(3) + t12;
t27 = t12 * t31;
t26 = t12 * t30;
t25 = t18 * t31;
t23 = t28 * t33;
t22 = t29 * t32;
t21 = -t6 * t12 - t18 * t24;
t20 = -t12 * t37 + t18 * t30;
t11 = t12 ^ 2;
t9 = -t18 * pkin(2) - pkin(3);
t8 = t16 * pkin(2) + pkin(6);
t4 = 0.2e1 * t15 * t26;
t2 = t6 * t31;
t1 = -0.2e1 * t34 * t12 * qJD(4);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t10; 0, 0, 0, 0, 0, t16 * t22, t18 * t22, t4, t1, t10, -t35, 0, t9 * t27 - t8 * t10 + t2 + (t29 * t36 - t25) * t32, -t20 * t32 + t26 * t9 + t35 * t8 + t38; 0, 0, 0, 0, 0, t16 * t23, t18 * t23, t4, t1, t10, -t35, 0, -pkin(3) * t27 - pkin(6) * t10 + t2 + (t28 * t36 + t25) * t33, -pkin(3) * t26 + pkin(6) * t35 + t20 * t33 + t38; 0, 0, 0, 0, 0, 0, 0, -t15 * t11 * t17, t34 * t11, 0, 0, 0, t21 * t15, t21 * t17;];
tauc_reg = t3;
