% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:40
% EndTime: 2019-12-31 17:34:40
% DurationCPUTime: 0.10s
% Computational Cost: add. (59->25), mult. (169->62), div. (0->0), fcn. (75->4), ass. (0->29)
t20 = qJD(3) ^ 2;
t29 = t20 / 0.2e1;
t16 = sin(qJ(4));
t28 = qJD(3) * t16;
t18 = cos(qJ(4));
t27 = qJD(3) * t18;
t26 = t18 * qJD(1);
t19 = cos(qJ(3));
t25 = t19 * qJD(2);
t24 = qJ(5) * qJD(3);
t23 = qJD(2) * qJD(3);
t22 = qJD(3) * qJD(4);
t17 = sin(qJ(3));
t7 = qJD(3) * pkin(6) + t17 * qJD(2);
t4 = -t16 * qJD(1) + t18 * t7;
t21 = qJD(2) ^ 2;
t15 = qJD(1) ^ 2 / 0.2e1;
t14 = qJD(4) ^ 2 / 0.2e1;
t13 = t18 * t22;
t12 = t16 * t22;
t11 = t18 ^ 2 * t29;
t10 = t16 ^ 2 * t29;
t9 = t16 * t20 * t18;
t8 = -qJD(3) * pkin(3) - t25;
t5 = -t25 + qJD(5) + (-pkin(4) * t18 - pkin(3)) * qJD(3);
t3 = -t16 * t7 - t26;
t2 = t18 * t24 + t4;
t1 = qJD(4) * pkin(4) - t26 + (-t7 - t24) * t16;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 + t21 / 0.2e1, 0, 0, 0, 0, 0, t29, t19 * t23, -t17 * t23, 0, t15 + (t17 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1) * t21, t10, t9, t12, t11, t13, t14, t3 * qJD(4) - t8 * t27, -t4 * qJD(4) + t8 * t28, (-t16 * t3 + t18 * t4) * qJD(3), t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1, t10, t9, t12, t11, t13, t14, t1 * qJD(4) - t5 * t27, -t2 * qJD(4) + t5 * t28, (-t1 * t16 + t18 * t2) * qJD(3), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
