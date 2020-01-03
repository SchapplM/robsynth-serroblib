% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:05
% EndTime: 2019-12-31 19:06:05
% DurationCPUTime: 0.12s
% Computational Cost: add. (237->33), mult. (349->87), div. (0->0), fcn. (149->6), ass. (0->35)
t27 = sin(qJ(4));
t24 = t27 ^ 2;
t40 = t24 / 0.2e1;
t29 = cos(qJ(4));
t25 = t29 ^ 2;
t39 = t25 / 0.2e1;
t31 = qJD(1) ^ 2;
t23 = t31 / 0.2e1;
t38 = cos(qJ(5));
t22 = -qJD(1) + qJD(3);
t37 = t22 * t27;
t36 = t22 * t29;
t16 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t33 = qJ(2) * qJD(1);
t13 = t28 * t16 + t30 * t33;
t35 = qJD(4) * t27;
t34 = qJD(4) * t29;
t8 = t22 * pkin(7) + t13;
t32 = pkin(8) * t22 + t8;
t12 = t30 * t16 - t28 * t33;
t26 = sin(qJ(5));
t21 = qJD(4) + qJD(5);
t20 = t22 ^ 2;
t19 = -pkin(1) * qJD(1) + qJD(2);
t11 = (t26 * t29 + t38 * t27) * t22;
t9 = t26 * t37 - t38 * t36;
t7 = -t22 * pkin(3) - t12;
t5 = (-pkin(4) * t29 - pkin(3)) * t22 - t12;
t4 = t32 * t29;
t3 = qJD(4) * pkin(4) - t32 * t27;
t2 = t26 * t3 + t38 * t4;
t1 = -t26 * t4 + t38 * t3;
t6 = [0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, -t19 * qJD(1), 0, t31 * qJ(2), qJ(2) ^ 2 * t23 + t19 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t20 / 0.2e1, t12 * t22, -t13 * t22, 0, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t20 * t40, t27 * t20 * t29, t22 * t35, t20 * t39, t22 * t34, qJD(4) ^ 2 / 0.2e1, -t8 * t35 - t7 * t36, -t8 * t34 + t7 * t37, (t24 + t25) * t8 * t22, t7 ^ 2 / 0.2e1 + (t39 + t40) * t8 ^ 2, t11 ^ 2 / 0.2e1, -t11 * t9, t11 * t21, t9 ^ 2 / 0.2e1, -t9 * t21, t21 ^ 2 / 0.2e1, t1 * t21 + t5 * t9, t5 * t11 - t2 * t21, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
