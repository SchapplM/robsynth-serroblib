% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:00
% EndTime: 2019-12-31 17:44:00
% DurationCPUTime: 0.12s
% Computational Cost: add. (127->36), mult. (328->81), div. (0->0), fcn. (168->6), ass. (0->32)
t29 = qJD(1) ^ 2;
t22 = t29 / 0.2e1;
t36 = pkin(1) * t29;
t24 = sin(pkin(7));
t18 = (pkin(1) * t24 + qJ(3)) * qJD(1);
t23 = sin(pkin(8));
t25 = cos(pkin(8));
t11 = t23 * qJD(2) + t25 * t18;
t35 = qJD(1) * t23;
t34 = qJD(1) * t25;
t33 = t23 * t29 * t25;
t26 = cos(pkin(7));
t32 = -pkin(1) * t26 - pkin(2);
t10 = t25 * qJD(2) - t23 * t18;
t9 = qJD(4) - t10;
t31 = qJ(4) * t23 - t32;
t28 = cos(qJ(5));
t27 = sin(qJ(5));
t20 = t25 ^ 2 * t22;
t19 = t23 ^ 2 * t22;
t17 = t32 * qJD(1) + qJD(3);
t14 = (t23 * t28 - t25 * t27) * qJD(1);
t12 = (-t23 * t27 - t25 * t28) * qJD(1);
t8 = t11 ^ 2 / 0.2e1;
t7 = qJD(3) + (-pkin(3) * t25 - t31) * qJD(1);
t6 = t11 * t34;
t5 = -pkin(6) * t34 + t11;
t4 = -pkin(6) * t35 + t9;
t3 = -qJD(3) + ((pkin(3) + pkin(4)) * t25 + t31) * qJD(1);
t2 = t27 * t4 + t28 * t5;
t1 = -t27 * t5 + t28 * t4;
t13 = [0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t26 * t36, -t24 * t36, 0, qJD(2) ^ 2 / 0.2e1 + (t24 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t29, t19, t33, 0, t20, 0, 0, -t17 * t34, t17 * t35, -t10 * t35 + t6, t8 + t10 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t19, 0, -t33, 0, 0, t20, -t7 * t34, t9 * t35 + t6, -t7 * t35, t8 + t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, t14 * t12, t14 * qJD(5), t12 ^ 2 / 0.2e1, t12 * qJD(5), qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t3 * t12, -t2 * qJD(5) + t3 * t14, -t1 * t14 + t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t13;
