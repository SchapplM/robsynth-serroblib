% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:01
% EndTime: 2019-12-31 17:51:01
% DurationCPUTime: 0.10s
% Computational Cost: add. (95->32), mult. (223->65), div. (0->0), fcn. (79->4), ass. (0->30)
t23 = qJD(1) ^ 2;
t18 = t23 / 0.2e1;
t34 = pkin(1) * t23;
t21 = sin(qJ(4));
t22 = cos(qJ(4));
t20 = cos(pkin(7));
t27 = -pkin(1) * t20 - pkin(2);
t8 = qJD(3) + (-pkin(6) + t27) * qJD(1);
t4 = t22 * qJD(2) + t21 * t8;
t19 = sin(pkin(7));
t25 = -pkin(1) * t19 - qJ(3);
t7 = qJD(5) + (pkin(4) * t21 - t25) * qJD(1);
t33 = qJD(1) * t7;
t10 = t25 * qJD(1);
t32 = t10 * qJD(1);
t31 = t10 ^ 2 / 0.2e1;
t30 = qJ(5) * qJD(1);
t29 = qJD(1) * qJD(4);
t28 = t22 * t23 * t21;
t26 = t21 * t29;
t3 = -t21 * qJD(2) + t22 * t8;
t17 = qJD(2) ^ 2 / 0.2e1;
t16 = qJD(4) ^ 2 / 0.2e1;
t14 = t22 * t29;
t13 = t22 ^ 2 * t18;
t12 = t21 ^ 2 * t18;
t9 = t27 * qJD(1) + qJD(3);
t2 = -t21 * t30 + t4;
t1 = qJD(4) * pkin(4) - t22 * t30 + t3;
t5 = [0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t20 * t34, -t19 * t34, 0, t17 + (t19 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t23, t18, 0, 0, 0, 0, 0, 0, t9 * qJD(1), -t32, t17 + t31 + t9 ^ 2 / 0.2e1, t13, -t28, t14, t12, -t26, t16, t3 * qJD(4) - t21 * t32, -t4 * qJD(4) - t22 * t32, (-t21 * t4 - t22 * t3) * qJD(1), t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t31, t13, -t28, t14, t12, -t26, t16, t1 * qJD(4) + t21 * t33, -t2 * qJD(4) + t22 * t33, (-t1 * t22 - t2 * t21) * qJD(1), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t5;
