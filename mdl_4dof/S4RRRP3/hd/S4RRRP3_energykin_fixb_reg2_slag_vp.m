% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRP3
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:14
% EndTime: 2019-12-31 17:14:14
% DurationCPUTime: 0.10s
% Computational Cost: add. (92->23), mult. (174->62), div. (0->0), fcn. (58->4), ass. (0->30)
t15 = sin(qJ(3));
t13 = t15 ^ 2;
t31 = t13 / 0.2e1;
t17 = cos(qJ(3));
t14 = t17 ^ 2;
t30 = t14 / 0.2e1;
t11 = qJD(1) + qJD(2);
t29 = t11 * t15;
t28 = t11 * t17;
t27 = pkin(1) * qJD(1);
t26 = qJD(3) * t15;
t25 = qJD(3) * t17;
t10 = t11 ^ 2;
t24 = t15 * t10 * t17;
t16 = sin(qJ(2));
t23 = t16 * t27;
t18 = cos(qJ(2));
t22 = t18 * t27;
t21 = t11 * t25;
t19 = qJD(1) ^ 2;
t12 = qJD(3) ^ 2 / 0.2e1;
t9 = t11 * t26;
t8 = t10 * t30;
t7 = t10 * t31;
t6 = -t11 * pkin(2) - t22;
t5 = t11 * pkin(6) + t23;
t3 = qJD(3) * qJ(4) + t17 * t5;
t2 = -qJD(3) * pkin(3) + t15 * t5 + qJD(4);
t1 = -t22 + (-pkin(3) * t17 - qJ(4) * t15 - pkin(2)) * t11;
t4 = [0, 0, 0, 0, 0, t19 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 / 0.2e1, t11 * t22, -t11 * t23, 0, (t16 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t19, t7, t24, t9, t8, t21, t12, -t5 * t26 - t6 * t28, -t5 * t25 + t6 * t29, (t13 + t14) * t5 * t11, t6 ^ 2 / 0.2e1 + (t30 + t31) * t5 ^ 2, t7, t9, -t24, t12, -t21, t8, -t2 * qJD(3) - t1 * t28, (t15 * t2 + t17 * t3) * t11, t3 * qJD(3) - t1 * t29, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t4;
