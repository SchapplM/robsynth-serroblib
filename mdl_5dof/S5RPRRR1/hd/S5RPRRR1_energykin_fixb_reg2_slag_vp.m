% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:09:50
% EndTime: 2019-12-05 18:09:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (93->31), mult. (306->87), div. (0->0), fcn. (190->6), ass. (0->31)
t27 = qJD(1) ^ 2;
t20 = t27 / 0.2e1;
t38 = cos(qJ(4));
t37 = cos(qJ(5));
t36 = t27 * qJ(2);
t24 = sin(qJ(4));
t26 = cos(qJ(3));
t34 = qJ(2) * qJD(1);
t7 = t24 * t26 * t34 - t38 * qJD(2);
t35 = t7 ^ 2 / 0.2e1;
t33 = qJ(2) * qJD(3);
t32 = qJD(1) * qJD(3);
t31 = qJ(2) ^ 2 * t20;
t25 = sin(qJ(3));
t30 = t25 * t34;
t29 = qJD(1) * t38;
t10 = t24 * t25 * qJD(1) - t38 * qJD(3);
t23 = sin(qJ(5));
t22 = t26 ^ 2;
t21 = t25 ^ 2;
t19 = qJD(2) ^ 2 / 0.2e1;
t15 = t26 * qJD(1) - qJD(4);
t14 = t21 * t31;
t12 = t24 * qJD(3) + t25 * t29;
t9 = t26 * qJ(2) * t29 + t24 * qJD(2);
t6 = qJD(5) + t10;
t5 = t23 * t30 + t37 * t9;
t4 = -t23 * t9 + t37 * t30;
t3 = t37 * t12 - t23 * t15;
t1 = t23 * t12 + t37 * t15;
t2 = [0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, -qJD(2) * qJD(1), 0, t36, t31 + t19, t21 * t20, t25 * t27 * t26, t25 * t32, t22 * t20, t26 * t32, qJD(3) ^ 2 / 0.2e1, (-qJD(2) * t26 - t25 * t33) * qJD(1), (qJD(2) * t25 - t26 * t33) * qJD(1), (t21 + t22) * t36, t22 * t31 + t14 + t19, t12 ^ 2 / 0.2e1, -t12 * t10, -t12 * t15, t10 ^ 2 / 0.2e1, t10 * t15, t15 ^ 2 / 0.2e1, t10 * t30 + t7 * t15, t12 * t30 + t9 * t15, -t9 * t10 + t7 * t12, t9 ^ 2 / 0.2e1 + t35 + t14, t3 ^ 2 / 0.2e1, -t3 * t1, t3 * t6, t1 ^ 2 / 0.2e1, -t1 * t6, t6 ^ 2 / 0.2e1, t7 * t1 + t4 * t6, t7 * t3 - t5 * t6, -t5 * t1 - t4 * t3, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t35;];
T_reg = t2;
