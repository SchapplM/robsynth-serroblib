% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:28
% EndTime: 2019-12-31 18:47:29
% DurationCPUTime: 0.11s
% Computational Cost: add. (162->30), mult. (249->69), div. (0->0), fcn. (80->4), ass. (0->33)
t23 = sin(qJ(4));
t21 = t23 ^ 2;
t36 = t21 / 0.2e1;
t25 = cos(qJ(4));
t22 = t25 ^ 2;
t35 = t22 / 0.2e1;
t27 = qJD(1) ^ 2;
t20 = t27 / 0.2e1;
t13 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t30 = qJ(2) * qJD(1);
t8 = t24 * t13 + t26 * t30;
t18 = -qJD(1) + qJD(3);
t34 = t18 * t23;
t33 = t18 * t25;
t32 = qJD(4) * t23;
t31 = qJD(4) * t25;
t17 = t18 ^ 2;
t29 = t23 * t17 * t25;
t28 = t18 * t31;
t7 = t26 * t13 - t24 * t30;
t19 = qJD(4) ^ 2 / 0.2e1;
t16 = -qJD(1) * pkin(1) + qJD(2);
t12 = t18 * t32;
t11 = t17 * t35;
t10 = t17 * t36;
t6 = t18 * pkin(7) + t8;
t5 = -t18 * pkin(3) - t7;
t3 = qJD(4) * qJ(5) + t25 * t6;
t2 = -qJD(4) * pkin(4) + t23 * t6 + qJD(5);
t1 = (-pkin(4) * t25 - qJ(5) * t23 - pkin(3)) * t18 - t7;
t4 = [0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, -t16 * qJD(1), 0, t27 * qJ(2), qJ(2) ^ 2 * t20 + t16 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t17 / 0.2e1, t7 * t18, -t8 * t18, 0, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t10, t29, t12, t11, t28, t19, -t6 * t32 - t5 * t33, -t6 * t31 + t5 * t34, (t21 + t22) * t6 * t18, t5 ^ 2 / 0.2e1 + (t35 + t36) * t6 ^ 2, t10, t12, -t29, t19, -t28, t11, -t2 * qJD(4) - t1 * t33, (t2 * t23 + t25 * t3) * t18, t3 * qJD(4) - t1 * t34, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t4;
