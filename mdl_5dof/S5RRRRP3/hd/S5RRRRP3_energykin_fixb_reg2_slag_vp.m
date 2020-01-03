% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:24
% EndTime: 2019-12-31 21:49:24
% DurationCPUTime: 0.13s
% Computational Cost: add. (208->29), mult. (300->75), div. (0->0), fcn. (128->6), ass. (0->36)
t22 = sin(qJ(4));
t20 = t22 ^ 2;
t40 = t20 / 0.2e1;
t25 = cos(qJ(4));
t21 = t25 ^ 2;
t39 = t21 / 0.2e1;
t18 = qJD(1) + qJD(2);
t27 = cos(qJ(2));
t36 = pkin(1) * qJD(1);
t31 = t27 * t36;
t10 = t18 * pkin(2) + t31;
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t24 = sin(qJ(2));
t32 = t24 * t36;
t8 = t23 * t10 + t26 * t32;
t17 = qJD(3) + t18;
t38 = t17 * t22;
t37 = t17 * t25;
t35 = qJD(4) * t22;
t34 = qJD(4) * t25;
t16 = t17 ^ 2;
t33 = t22 * t16 * t25;
t30 = t17 * t34;
t7 = t26 * t10 - t23 * t32;
t28 = qJD(1) ^ 2;
t19 = qJD(4) ^ 2 / 0.2e1;
t13 = t17 * t35;
t12 = t16 * t39;
t11 = t16 * t40;
t6 = t17 * pkin(8) + t8;
t5 = -t17 * pkin(3) - t7;
t3 = qJD(4) * qJ(5) + t25 * t6;
t2 = -qJD(4) * pkin(4) + t22 * t6 + qJD(5);
t1 = (-pkin(4) * t25 - qJ(5) * t22 - pkin(3)) * t17 - t7;
t4 = [0, 0, 0, 0, 0, t28 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 ^ 2 / 0.2e1, t18 * t31, -t18 * t32, 0, (t24 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t28, 0, 0, 0, 0, 0, t16 / 0.2e1, t7 * t17, -t8 * t17, 0, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t11, t33, t13, t12, t30, t19, -t6 * t35 - t5 * t37, -t6 * t34 + t5 * t38, (t20 + t21) * t6 * t17, t5 ^ 2 / 0.2e1 + (t39 + t40) * t6 ^ 2, t11, t13, -t33, t19, -t30, t12, -t2 * qJD(4) - t1 * t37, (t2 * t22 + t25 * t3) * t17, t3 * qJD(4) - t1 * t38, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t4;
