% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:03
% EndTime: 2019-12-05 15:31:03
% DurationCPUTime: 0.15s
% Computational Cost: add. (119->33), mult. (317->73), div. (0->0), fcn. (169->4), ass. (0->35)
t27 = qJD(2) ^ 2;
t40 = t27 / 0.2e1;
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t35 = qJ(3) * qJD(2);
t13 = t23 * qJD(1) + t24 * t35;
t25 = sin(qJ(4));
t26 = cos(qJ(4));
t8 = qJD(3) + (-pkin(3) * t24 - pkin(6) * t23 - pkin(2)) * qJD(2);
t4 = t26 * t13 + t25 * t8;
t39 = t23 ^ 2 * t27;
t38 = qJD(2) * t23;
t37 = t24 * qJD(2);
t20 = t24 * qJD(1);
t11 = t23 * t35 - t20;
t36 = t11 ^ 2 / 0.2e1;
t34 = t25 * t38;
t33 = t26 * t38;
t32 = t11 * t38;
t31 = t39 / 0.2e1;
t30 = qJ(5) * t38;
t3 = -t25 * t13 + t26 * t8;
t29 = t26 * t25 * t39;
t17 = -qJD(4) + t37;
t28 = t17 * t34;
t22 = qJD(1) ^ 2 / 0.2e1;
t19 = -qJD(2) * pkin(2) + qJD(3);
t16 = t26 ^ 2 * t31;
t15 = t25 ^ 2 * t31;
t14 = t17 ^ 2 / 0.2e1;
t10 = t17 * t33;
t5 = qJD(5) - t20 + (pkin(4) * t25 + qJ(3)) * t38;
t2 = -t25 * t30 + t4;
t1 = -t17 * pkin(4) - t26 * t30 + t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, t40, 0, 0, 0, t22, t31, t23 * t27 * t24, 0, t24 ^ 2 * t40, 0, 0, -t19 * t37, t19 * t38, (t11 * t23 + t13 * t24) * qJD(2), t13 ^ 2 / 0.2e1 + t36 + t19 ^ 2 / 0.2e1, t16, -t29, -t10, t15, t28, t14, -t3 * t17 + t25 * t32, t4 * t17 + t26 * t32, (-t25 * t4 - t26 * t3) * t38, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t36, t16, -t29, -t10, t15, t28, t14, -t1 * t17 + t5 * t34, t2 * t17 + t5 * t33, (-t1 * t26 - t2 * t25) * t38, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
