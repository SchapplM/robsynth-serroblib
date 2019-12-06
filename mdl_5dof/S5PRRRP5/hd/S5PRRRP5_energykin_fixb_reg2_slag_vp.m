% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:16
% EndTime: 2019-12-05 16:49:16
% DurationCPUTime: 0.14s
% Computational Cost: add. (167->37), mult. (425->91), div. (0->0), fcn. (251->6), ass. (0->38)
t33 = qJD(2) ^ 2;
t43 = t33 / 0.2e1;
t29 = sin(qJ(3));
t30 = sin(qJ(2));
t21 = qJD(2) * pkin(6) + t30 * qJD(1);
t35 = pkin(7) * qJD(2) + t21;
t10 = qJD(3) * pkin(3) - t35 * t29;
t31 = cos(qJ(3));
t13 = t35 * t31;
t28 = sin(qJ(4));
t42 = cos(qJ(4));
t4 = t28 * t10 + t42 * t13;
t41 = qJD(2) * t29;
t40 = qJD(2) * t31;
t39 = qJD(3) * t21;
t32 = cos(qJ(2));
t38 = t32 * qJD(1);
t37 = qJD(1) * qJD(2);
t36 = qJD(2) * qJD(3);
t3 = t42 * t10 - t28 * t13;
t19 = -t38 + (-pkin(3) * t31 - pkin(2)) * qJD(2);
t34 = qJD(1) ^ 2;
t27 = t31 ^ 2;
t26 = t29 ^ 2;
t25 = qJD(3) + qJD(4);
t24 = t25 ^ 2 / 0.2e1;
t22 = -qJD(2) * pkin(2) - t38;
t18 = (t28 * t31 + t42 * t29) * qJD(2);
t16 = t28 * t41 - t42 * t40;
t15 = t18 ^ 2 / 0.2e1;
t14 = t16 ^ 2 / 0.2e1;
t12 = t18 * t25;
t11 = t16 * t25;
t6 = t18 * t16;
t5 = t16 * pkin(4) + qJD(5) + t19;
t2 = -t16 * qJ(5) + t4;
t1 = t25 * pkin(4) - t18 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t34 / 0.2e1, 0, 0, 0, 0, 0, t43, t32 * t37, -t30 * t37, 0, (t30 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1) * t34, t26 * t43, t29 * t33 * t31, t29 * t36, t27 * t43, t31 * t36, qJD(3) ^ 2 / 0.2e1, -t22 * t40 - t29 * t39, t22 * t41 - t31 * t39, (t26 + t27) * t21 * qJD(2), t22 ^ 2 / 0.2e1 + (t27 / 0.2e1 + t26 / 0.2e1) * t21 ^ 2, t15, -t6, t12, t14, -t11, t24, t19 * t16 + t3 * t25, t19 * t18 - t4 * t25, -t4 * t16 - t3 * t18, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t15, -t6, t12, t14, -t11, t24, t1 * t25 + t5 * t16, t5 * t18 - t2 * t25, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
