% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:25:31
% EndTime: 2019-12-05 18:25:31
% DurationCPUTime: 0.16s
% Computational Cost: add. (199->37), mult. (501->91), div. (0->0), fcn. (288->6), ass. (0->39)
t36 = qJD(1) ^ 2;
t47 = t36 / 0.2e1;
t46 = pkin(2) + pkin(1);
t45 = cos(qJ(4));
t44 = cos(qJ(5));
t35 = cos(qJ(2));
t43 = t35 ^ 2 * t36;
t42 = pkin(3) + qJ(3);
t39 = qJD(1) * t35;
t19 = t42 * t39;
t33 = sin(qJ(4));
t34 = sin(qJ(2));
t40 = qJD(1) * t34;
t37 = t46 * qJD(2) - t42 * t40;
t6 = t45 * t19 + t33 * t37;
t4 = t33 * t19 - t45 * t37;
t41 = t4 ^ 2 / 0.2e1;
t38 = qJD(1) * qJD(2);
t26 = t43 / 0.2e1;
t15 = t33 * t40 - t45 * t39;
t18 = -t46 * t39 + qJD(3);
t32 = sin(qJ(5));
t30 = qJD(2) ^ 2 / 0.2e1;
t29 = qJD(2) + qJD(4);
t28 = t35 * t38;
t27 = t34 * t38;
t25 = t34 ^ 2 * t47;
t22 = t34 * t36 * t35;
t21 = -pkin(1) * t39 + qJD(3);
t20 = qJD(2) * pkin(1) - qJ(3) * t40;
t17 = (t33 * t35 + t45 * t34) * qJD(1);
t13 = qJD(5) + t15;
t10 = t44 * t17 + t32 * t29;
t8 = t32 * t17 - t44 * t29;
t7 = -t17 * pkin(4) + t18;
t3 = t29 * pkin(4) + t6;
t2 = t44 * t3 + t32 * t7;
t1 = -t32 * t3 + t44 * t7;
t5 = [0, 0, 0, 0, 0, t47, 0, 0, 0, 0, t25, t22, t27, t26, t28, t30, 0, 0, 0, 0, t25, t22, t27, t26, t28, t30, t20 * qJD(2) - t21 * t39, (-qJ(3) * qJD(2) * t35 + t21 * t34) * qJD(1), qJ(3) * t43 - t20 * t40, qJ(3) ^ 2 * t26 + t20 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t15, t17 * t29, t15 ^ 2 / 0.2e1, -t15 * t29, t29 ^ 2 / 0.2e1, t18 * t15 - t4 * t29, t18 * t17 - t6 * t29, -t6 * t15 + t4 * t17, t6 ^ 2 / 0.2e1 + t41 + t18 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t13, t8 ^ 2 / 0.2e1, -t8 * t13, t13 ^ 2 / 0.2e1, t1 * t13 + t4 * t8, t4 * t10 - t2 * t13, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t41;];
T_reg = t5;
