% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:05
% EndTime: 2019-12-31 17:17:05
% DurationCPUTime: 0.11s
% Computational Cost: add. (126->30), mult. (363->77), div. (0->0), fcn. (196->4), ass. (0->33)
t24 = qJD(1) ^ 2;
t37 = t24 / 0.2e1;
t36 = -pkin(6) - pkin(5);
t21 = sin(qJ(3));
t22 = sin(qJ(2));
t23 = cos(qJ(2));
t33 = cos(qJ(3));
t11 = (t21 * t23 + t33 * t22) * qJD(1);
t29 = qJD(1) * t23;
t30 = qJD(1) * t22;
t9 = t21 * t30 - t33 * t29;
t35 = t11 * t9;
t18 = qJD(2) + qJD(3);
t34 = t18 * t9;
t13 = qJD(2) * pkin(2) + t36 * t30;
t14 = t36 * t29;
t5 = t21 * t13 - t33 * t14;
t32 = t23 * t24;
t31 = t9 ^ 2 / 0.2e1;
t28 = qJD(1) * qJD(2);
t27 = t22 * t28;
t26 = t23 * t28;
t15 = (-pkin(2) * t23 - pkin(1)) * qJD(1);
t4 = t33 * t13 + t21 * t14;
t20 = t23 ^ 2;
t19 = t22 ^ 2;
t17 = t18 ^ 2 / 0.2e1;
t8 = t11 ^ 2 / 0.2e1;
t6 = t11 * t18;
t3 = t18 * qJ(4) + t5;
t2 = -t18 * pkin(3) + qJD(4) - t4;
t1 = t9 * pkin(3) - t11 * qJ(4) + t15;
t7 = [0, 0, 0, 0, 0, t37, 0, 0, 0, 0, t19 * t37, t22 * t32, t27, t20 * t37, t26, qJD(2) ^ 2 / 0.2e1, pkin(1) * t32 - pkin(5) * t27, -t24 * pkin(1) * t22 - pkin(5) * t26, (t19 + t20) * t24 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t20 / 0.2e1 + t19 / 0.2e1) * pkin(5) ^ 2) * t24, t8, -t35, t6, t31, -t34, t17, t15 * t9 + t4 * t18, t15 * t11 - t5 * t18, -t4 * t11 - t5 * t9, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t8, t6, t35, t17, t34, t31, t1 * t9 - t2 * t18, t2 * t11 - t3 * t9, -t1 * t11 + t3 * t18, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t7;
