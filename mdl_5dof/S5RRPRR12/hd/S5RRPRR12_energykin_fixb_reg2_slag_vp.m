% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:36
% EndTime: 2019-12-31 20:30:36
% DurationCPUTime: 0.19s
% Computational Cost: add. (272->50), mult. (622->116), div. (0->0), fcn. (333->6), ass. (0->41)
t40 = sin(qJ(4));
t41 = sin(qJ(2));
t42 = cos(qJ(4));
t43 = cos(qJ(2));
t16 = (t40 * t41 + t42 * t43) * qJD(1);
t44 = qJD(1) ^ 2;
t55 = t44 / 0.2e1;
t54 = cos(qJ(5));
t53 = t43 * t44;
t52 = qJD(1) * t41;
t50 = pkin(6) * t52 + qJD(3);
t13 = -pkin(7) * t52 + (-pkin(2) - pkin(3)) * qJD(2) + t50;
t51 = qJD(1) * t43;
t22 = pkin(6) * t51 + qJD(2) * qJ(3);
t19 = -pkin(7) * t51 + t22;
t7 = t40 * t13 + t42 * t19;
t49 = qJD(1) * qJD(2);
t48 = t41 * t53;
t20 = -qJD(1) * pkin(1) - pkin(2) * t51 - qJ(3) * t52;
t26 = t41 * t49;
t47 = t43 * t49;
t12 = pkin(3) * t51 - t20;
t6 = t42 * t13 - t40 * t19;
t39 = sin(qJ(5));
t37 = t43 ^ 2;
t36 = t41 ^ 2;
t34 = qJD(2) ^ 2 / 0.2e1;
t32 = qJD(2) - qJD(4);
t25 = t37 * t55;
t24 = t36 * t55;
t21 = -qJD(2) * pkin(2) + t50;
t18 = (-t40 * t43 + t41 * t42) * qJD(1);
t15 = qJD(5) + t16;
t10 = t54 * t18 - t39 * t32;
t8 = t39 * t18 + t54 * t32;
t5 = -t32 * pkin(8) + t7;
t4 = t32 * pkin(4) - t6;
t3 = t16 * pkin(4) - t18 * pkin(8) + t12;
t2 = t39 * t3 + t54 * t5;
t1 = t54 * t3 - t39 * t5;
t9 = [0, 0, 0, 0, 0, t55, 0, 0, 0, 0, t24, t48, t26, t25, t47, t34, pkin(1) * t53 - pkin(6) * t26, -t44 * pkin(1) * t41 - pkin(6) * t47, (t36 + t37) * t44 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t37 / 0.2e1 + t36 / 0.2e1) * pkin(6) ^ 2) * t44, t24, t26, -t48, t34, -t47, t25, -t21 * qJD(2) - t20 * t51, (t21 * t41 + t22 * t43) * qJD(1), t22 * qJD(2) - t20 * t52, t22 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, -t18 * t32, t16 ^ 2 / 0.2e1, t16 * t32, t32 ^ 2 / 0.2e1, t12 * t16 - t6 * t32, t12 * t18 + t7 * t32, -t7 * t16 - t6 * t18, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t15, t8 ^ 2 / 0.2e1, -t8 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t4 * t8, t4 * t10 - t2 * t15, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg = t9;
