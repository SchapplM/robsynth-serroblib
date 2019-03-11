% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:37:45
% EndTime: 2019-03-08 18:37:45
% DurationCPUTime: 0.16s
% Computational Cost: add. (309->42), mult. (757->118), div. (0->0), fcn. (564->8), ass. (0->36)
t36 = qJD(1) ^ 2;
t45 = t36 / 0.2e1;
t44 = cos(qJ(4));
t43 = cos(qJ(5));
t34 = cos(qJ(2));
t42 = t34 * t36;
t27 = qJD(2) + qJD(3);
t33 = cos(qJ(3));
t41 = pkin(2) * qJD(2);
t38 = t33 * t41;
t21 = t27 * pkin(3) + t38;
t30 = sin(qJ(4));
t31 = sin(qJ(3));
t39 = t31 * t41;
t15 = t30 * t21 + t44 * t39;
t22 = (pkin(2) * t34 + pkin(1)) * qJD(1);
t40 = qJD(1) * qJD(2);
t32 = sin(qJ(2));
t17 = (t31 * t32 - t33 * t34) * qJD(1);
t18 = (-t31 * t34 - t32 * t33) * qJD(1);
t8 = -t44 * t17 + t30 * t18;
t13 = -t17 * pkin(3) + t22;
t14 = t44 * t21 - t30 * t39;
t35 = qJD(2) ^ 2;
t29 = sin(qJ(5));
t26 = qJD(4) + t27;
t12 = t26 * pkin(6) + t15;
t11 = -t26 * pkin(4) - t14;
t10 = t30 * t17 + t44 * t18;
t7 = qJD(5) + t8;
t6 = t43 * t10 + t29 * t26;
t4 = t29 * t10 - t43 * t26;
t3 = t8 * pkin(4) - t10 * pkin(6) + t13;
t2 = t43 * t12 + t29 * t3;
t1 = -t29 * t12 + t43 * t3;
t5 = [0, 0, 0, 0, 0, t45, 0, 0, 0, 0, t32 ^ 2 * t45, t32 * t42, -t32 * t40, t34 ^ 2 * t45, -t34 * t40, t35 / 0.2e1, pkin(1) * t42, -t36 * pkin(1) * t32, 0, pkin(1) ^ 2 * t45, t18 ^ 2 / 0.2e1, t18 * t17, t18 * t27, t17 ^ 2 / 0.2e1, t17 * t27, t27 ^ 2 / 0.2e1, -t22 * t17 + t27 * t38, t22 * t18 - t27 * t39 (t17 * t31 - t18 * t33) * t41, t22 ^ 2 / 0.2e1 + (t31 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t35, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t26, t8 ^ 2 / 0.2e1, -t8 * t26, t26 ^ 2 / 0.2e1, t13 * t8 + t14 * t26, t13 * t10 - t15 * t26, -t14 * t10 - t15 * t8, t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t6 ^ 2 / 0.2e1, -t6 * t4, t6 * t7, t4 ^ 2 / 0.2e1, -t4 * t7, t7 ^ 2 / 0.2e1, t1 * t7 + t11 * t4, t11 * t6 - t2 * t7, -t1 * t6 - t2 * t4, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1;];
T_reg  = t5;
