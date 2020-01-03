% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:43
% EndTime: 2019-12-31 17:06:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (172->36), mult. (490->94), div. (0->0), fcn. (293->6), ass. (0->34)
t29 = qJD(1) ^ 2;
t40 = t29 / 0.2e1;
t39 = cos(qJ(4));
t28 = cos(qJ(2));
t38 = t28 * t29;
t37 = pkin(5) + qJ(3);
t27 = sin(qJ(2));
t35 = qJD(1) * t27;
t17 = qJD(2) * pkin(2) - t37 * t35;
t34 = qJD(1) * t28;
t18 = t37 * t34;
t25 = sin(pkin(7));
t36 = cos(pkin(7));
t7 = t25 * t17 + t36 * t18;
t33 = qJD(1) * qJD(2);
t32 = t27 * t33;
t31 = t28 * t33;
t13 = t25 * t35 - t36 * t34;
t6 = t36 * t17 - t25 * t18;
t19 = qJD(3) + (-pkin(2) * t28 - pkin(1)) * qJD(1);
t26 = sin(qJ(4));
t24 = t28 ^ 2;
t23 = t27 ^ 2;
t22 = qJD(2) ^ 2 / 0.2e1;
t15 = (t25 * t28 + t36 * t27) * qJD(1);
t12 = qJD(4) + t13;
t10 = t26 * qJD(2) + t39 * t15;
t8 = -t39 * qJD(2) + t26 * t15;
t5 = qJD(2) * pkin(6) + t7;
t4 = -qJD(2) * pkin(3) - t6;
t3 = t13 * pkin(3) - t15 * pkin(6) + t19;
t2 = t26 * t3 + t39 * t5;
t1 = -t26 * t5 + t39 * t3;
t9 = [0, 0, 0, 0, 0, t40, 0, 0, 0, 0, t23 * t40, t27 * t38, t32, t24 * t40, t31, t22, pkin(1) * t38 - pkin(5) * t32, -t29 * pkin(1) * t27 - pkin(5) * t31, (t23 + t24) * t29 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t24 / 0.2e1 + t23 / 0.2e1) * pkin(5) ^ 2) * t29, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * qJD(2), t13 ^ 2 / 0.2e1, -t13 * qJD(2), t22, t6 * qJD(2) + t19 * t13, -t7 * qJD(2) + t19 * t15, -t7 * t13 - t6 * t15, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t12, t8 ^ 2 / 0.2e1, -t8 * t12, t12 ^ 2 / 0.2e1, t1 * t12 + t4 * t8, t4 * t10 - t2 * t12, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg = t9;
