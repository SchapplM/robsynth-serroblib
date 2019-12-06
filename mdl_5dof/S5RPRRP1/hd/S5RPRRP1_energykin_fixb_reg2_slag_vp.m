% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP1
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:56
% EndTime: 2019-12-05 17:59:56
% DurationCPUTime: 0.14s
% Computational Cost: add. (191->34), mult. (416->82), div. (0->0), fcn. (208->4), ass. (0->34)
t36 = qJD(1) ^ 2;
t28 = t36 / 0.2e1;
t35 = cos(qJ(3));
t21 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t37 = -pkin(7) * qJD(1) + t21;
t10 = qJD(3) * pkin(3) + t37 * t35;
t33 = sin(qJ(3));
t13 = t37 * t33;
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t4 = t32 * t10 + t34 * t13;
t40 = qJD(1) * t33;
t19 = pkin(3) * t40 + qJD(1) * qJ(2);
t41 = t36 * qJ(2);
t39 = qJD(3) * t21;
t38 = qJD(1) * qJD(3);
t3 = t34 * t10 - t32 * t13;
t31 = t35 ^ 2;
t30 = t33 ^ 2;
t27 = qJD(3) + qJD(4);
t26 = qJ(2) ^ 2 * t28;
t25 = -pkin(1) * qJD(1) + qJD(2);
t23 = t27 ^ 2 / 0.2e1;
t18 = t34 * t35 * qJD(1) - t32 * t40;
t16 = (-t32 * t35 - t33 * t34) * qJD(1);
t15 = t18 ^ 2 / 0.2e1;
t14 = t16 ^ 2 / 0.2e1;
t12 = t18 * t27;
t11 = t16 * t27;
t6 = -t16 * pkin(4) + qJD(5) + t19;
t5 = t18 * t16;
t2 = t16 * qJ(5) + t4;
t1 = t27 * pkin(4) - t18 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, t28, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, t25 * qJD(1), t41, t26 + t25 ^ 2 / 0.2e1, t31 * t28, -t35 * t36 * t33, t35 * t38, t30 * t28, -t33 * t38, qJD(3) ^ 2 / 0.2e1, t33 * t41 + t35 * t39, -t33 * t39 + t35 * t41, (-t30 - t31) * t21 * qJD(1), t26 + (t30 / 0.2e1 + t31 / 0.2e1) * t21 ^ 2, t15, t5, t12, t14, t11, t23, -t19 * t16 + t3 * t27, t19 * t18 - t4 * t27, t4 * t16 - t3 * t18, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t15, t5, t12, t14, t11, t23, t1 * t27 - t6 * t16, t6 * t18 - t2 * t27, -t1 * t18 + t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg = t7;
