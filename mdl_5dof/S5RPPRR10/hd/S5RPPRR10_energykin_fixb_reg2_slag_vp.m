% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:16
% EndTime: 2019-12-31 18:04:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (217->42), mult. (589->97), div. (0->0), fcn. (365->6), ass. (0->38)
t41 = qJD(1) ^ 2;
t49 = t41 / 0.2e1;
t48 = cos(qJ(5));
t35 = sin(pkin(8));
t46 = qJD(1) * t35;
t22 = qJ(2) * t46 + qJD(3);
t20 = -pkin(6) * t46 + t22;
t36 = cos(pkin(8));
t45 = qJD(1) * t36;
t21 = (-pkin(6) + qJ(2)) * t45;
t39 = sin(qJ(4));
t40 = cos(qJ(4));
t10 = t39 * t20 + t40 * t21;
t47 = qJ(2) * t41;
t31 = -qJD(1) * pkin(1) + qJD(2);
t44 = t35 * t41 * t36;
t43 = qJ(2) ^ 2 * t49;
t9 = t40 * t20 - t39 * t21;
t15 = -pkin(2) * t45 - qJ(3) * t46 + t31;
t11 = pkin(3) * t45 - t15;
t38 = sin(qJ(5));
t34 = qJD(4) + qJD(5);
t33 = t36 ^ 2;
t32 = t35 ^ 2;
t28 = t33 * t47;
t25 = t33 * t49;
t24 = t32 * t49;
t23 = t33 * t43;
t19 = (t35 * t40 - t36 * t39) * qJD(1);
t17 = (-t35 * t39 - t36 * t40) * qJD(1);
t8 = t38 * t17 + t48 * t19;
t6 = -t48 * t17 + t38 * t19;
t5 = -t17 * pkin(4) + t11;
t4 = t17 * pkin(7) + t10;
t3 = qJD(4) * pkin(4) - t19 * pkin(7) + t9;
t2 = t38 * t3 + t48 * t4;
t1 = t48 * t3 - t38 * t4;
t7 = [0, 0, 0, 0, 0, t49, 0, 0, 0, 0, t24, t44, 0, t25, 0, 0, -t31 * t45, t31 * t46, t32 * t47 + t28, t23 + t32 * t43 + t31 ^ 2 / 0.2e1, t24, 0, -t44, 0, 0, t25, -t15 * t45, t22 * t46 + t28, -t15 * t46, t23 + t15 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, t19 * t17, t19 * qJD(4), t17 ^ 2 / 0.2e1, t17 * qJD(4), qJD(4) ^ 2 / 0.2e1, t9 * qJD(4) - t11 * t17, -t10 * qJD(4) + t11 * t19, t10 * t17 - t9 * t19, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t8 ^ 2 / 0.2e1, -t8 * t6, t8 * t34, t6 ^ 2 / 0.2e1, -t6 * t34, t34 ^ 2 / 0.2e1, t1 * t34 + t5 * t6, -t2 * t34 + t5 * t8, -t1 * t8 - t2 * t6, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
