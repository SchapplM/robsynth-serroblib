% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:39:15
% EndTime: 2020-01-03 11:39:15
% DurationCPUTime: 0.16s
% Computational Cost: add. (282->45), mult. (707->114), div. (0->0), fcn. (439->8), ass. (0->38)
t38 = qJD(1) ^ 2;
t31 = t38 / 0.2e1;
t47 = pkin(1) * t38;
t46 = cos(qJ(5));
t33 = sin(pkin(8));
t24 = (pkin(1) * t33 + pkin(6)) * qJD(1);
t37 = cos(qJ(3));
t28 = t37 * qJD(2);
t36 = sin(qJ(3));
t42 = qJ(4) * qJD(1);
t14 = qJD(3) * pkin(3) + t28 + (-t24 - t42) * t36;
t18 = t36 * qJD(2) + t37 * t24;
t15 = t37 * t42 + t18;
t32 = sin(pkin(9));
t45 = cos(pkin(9));
t6 = t32 * t14 + t45 * t15;
t44 = qJD(1) * t36;
t43 = qJD(1) * t37;
t41 = qJD(1) * qJD(3);
t34 = cos(pkin(8));
t40 = -pkin(1) * t34 - pkin(2);
t5 = t45 * t14 - t32 * t15;
t19 = qJD(4) + (-pkin(3) * t37 + t40) * qJD(1);
t35 = sin(qJ(5));
t30 = qJD(3) ^ 2 / 0.2e1;
t29 = qJD(3) + qJD(5);
t25 = t40 * qJD(1);
t22 = (t32 * t37 + t45 * t36) * qJD(1);
t20 = t32 * t44 - t45 * t43;
t17 = -t36 * t24 + t28;
t10 = t20 * pkin(4) + t19;
t9 = -t35 * t20 + t46 * t22;
t7 = t46 * t20 + t35 * t22;
t4 = -t20 * pkin(7) + t6;
t3 = qJD(3) * pkin(4) - t22 * pkin(7) + t5;
t2 = t35 * t3 + t46 * t4;
t1 = t46 * t3 - t35 * t4;
t8 = [0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t34 * t47, -t33 * t47, 0, qJD(2) ^ 2 / 0.2e1 + (t33 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t38, t36 ^ 2 * t31, t36 * t38 * t37, t36 * t41, t37 ^ 2 * t31, t37 * t41, t30, t17 * qJD(3) - t25 * t43, -t18 * qJD(3) + t25 * t44, (-t17 * t36 + t18 * t37) * qJD(1), t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * qJD(3), t20 ^ 2 / 0.2e1, -t20 * qJD(3), t30, t5 * qJD(3) + t19 * t20, -t6 * qJD(3) + t19 * t22, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t29, t7 ^ 2 / 0.2e1, -t7 * t29, t29 ^ 2 / 0.2e1, t1 * t29 + t10 * t7, t10 * t9 - t2 * t29, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
