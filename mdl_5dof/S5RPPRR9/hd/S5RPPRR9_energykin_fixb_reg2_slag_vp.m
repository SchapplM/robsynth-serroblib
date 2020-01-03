% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:47
% EndTime: 2019-12-31 18:02:47
% DurationCPUTime: 0.12s
% Computational Cost: add. (184->37), mult. (353->89), div. (0->0), fcn. (155->6), ass. (0->29)
t31 = qJD(1) ^ 2;
t24 = t31 / 0.2e1;
t17 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t25 = sin(pkin(8));
t26 = cos(pkin(8));
t33 = qJ(2) * qJD(1);
t12 = t25 * t17 + t26 * t33;
t10 = -qJD(1) * pkin(6) + t12;
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t7 = t28 * qJD(3) + t30 * t10;
t36 = sin(qJ(5));
t35 = qJD(1) * t28;
t34 = t30 * qJD(1);
t32 = qJD(1) * qJD(4);
t11 = t26 * t17 - t25 * t33;
t9 = qJD(1) * pkin(3) - t11;
t6 = t30 * qJD(3) - t28 * t10;
t29 = cos(qJ(5));
t21 = -pkin(1) * qJD(1) + qJD(2);
t18 = qJD(5) + t34;
t14 = -t36 * qJD(4) + t29 * t35;
t13 = t29 * qJD(4) + t36 * t35;
t5 = (pkin(4) * t30 + pkin(7) * t28) * qJD(1) + t9;
t4 = qJD(4) * pkin(7) + t7;
t3 = -qJD(4) * pkin(4) - t6;
t2 = t29 * t4 + t36 * t5;
t1 = t29 * t5 - t36 * t4;
t8 = [0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, -t21 * qJD(1), 0, t31 * qJ(2), qJ(2) ^ 2 * t24 + t21 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t24, -t11 * qJD(1), t12 * qJD(1), 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t28 ^ 2 * t24, t28 * t31 * t30, -t28 * t32, t30 ^ 2 * t24, -t30 * t32, qJD(4) ^ 2 / 0.2e1, t6 * qJD(4) + t9 * t34, -t7 * qJD(4) - t9 * t35, (t28 * t6 - t30 * t7) * qJD(1), t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t13, -t14 * t18, t13 ^ 2 / 0.2e1, t13 * t18, t18 ^ 2 / 0.2e1, t1 * t18 - t3 * t13, -t3 * t14 - t2 * t18, t1 * t14 + t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
