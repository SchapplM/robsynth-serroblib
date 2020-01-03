% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP13
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP13_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:36
% EndTime: 2019-12-31 18:59:36
% DurationCPUTime: 0.12s
% Computational Cost: add. (179->37), mult. (376->82), div. (0->0), fcn. (174->4), ass. (0->33)
t29 = qJD(1) ^ 2;
t23 = t29 / 0.2e1;
t18 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t27 = sin(qJ(3));
t10 = qJD(3) * pkin(7) + t27 * t18;
t26 = sin(qJ(4));
t37 = cos(qJ(4));
t28 = cos(qJ(3));
t9 = (pkin(3) * t27 - pkin(7) * t28 + qJ(2)) * qJD(1);
t5 = t37 * t10 + t26 * t9;
t33 = qJD(1) * t28;
t13 = -t37 * qJD(3) + t26 * t33;
t15 = t26 * qJD(3) + t37 * t33;
t36 = t15 * t13;
t19 = t27 * qJD(1) + qJD(4);
t35 = t19 * t13;
t34 = t29 * qJ(2);
t32 = qJD(3) * t18;
t31 = t13 ^ 2 / 0.2e1;
t30 = qJD(1) * qJD(3);
t11 = -qJD(3) * pkin(3) - t28 * t18;
t4 = -t26 * t10 + t37 * t9;
t25 = t28 ^ 2;
t24 = t27 ^ 2;
t21 = qJ(2) ^ 2 * t23;
t20 = -pkin(1) * qJD(1) + qJD(2);
t16 = t19 ^ 2 / 0.2e1;
t12 = t15 ^ 2 / 0.2e1;
t6 = t15 * t19;
t3 = t13 * pkin(4) - t15 * qJ(5) + t11;
t2 = t19 * qJ(5) + t5;
t1 = -t19 * pkin(4) + qJD(5) - t4;
t7 = [0, 0, 0, 0, 0, t23, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, t20 * qJD(1), t34, t21 + t20 ^ 2 / 0.2e1, t25 * t23, -t28 * t29 * t27, t28 * t30, t24 * t23, -t27 * t30, qJD(3) ^ 2 / 0.2e1, t27 * t34 + t28 * t32, -t27 * t32 + t28 * t34, (-t24 - t25) * t18 * qJD(1), t21 + (t24 / 0.2e1 + t25 / 0.2e1) * t18 ^ 2, t12, -t36, t6, t31, -t35, t16, t11 * t13 + t4 * t19, t11 * t15 - t5 * t19, -t5 * t13 - t4 * t15, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t12, t6, t36, t16, t35, t31, -t1 * t19 + t3 * t13, t1 * t15 - t2 * t13, -t3 * t15 + t2 * t19, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t7;
