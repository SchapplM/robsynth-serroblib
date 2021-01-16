% Calculate minimal parameter regressor of fixed base kinetic energy for
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
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:45
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:44:13
% EndTime: 2021-01-15 11:44:13
% DurationCPUTime: 0.21s
% Computational Cost: add. (158->36), mult. (391->84), div. (0->0), fcn. (241->8), ass. (0->33)
t48 = qJD(1) ^ 2;
t57 = t48 / 0.2e1;
t42 = sin(pkin(8));
t35 = (pkin(1) * t42 + pkin(6)) * qJD(1);
t47 = cos(qJ(3));
t39 = t47 * qJD(2);
t45 = sin(qJ(3));
t52 = qJ(4) * qJD(1);
t28 = qJD(3) * pkin(3) + t39 + (-t35 - t52) * t45;
t56 = t45 * qJD(2) + t47 * t35;
t29 = t47 * t52 + t56;
t41 = sin(pkin(9));
t55 = cos(pkin(9));
t21 = t41 * t28 + t55 * t29;
t54 = qJD(1) * t45;
t53 = qJD(1) * t47;
t51 = qJD(1) * qJD(3);
t43 = cos(pkin(8));
t50 = -pkin(1) * t43 - pkin(2);
t20 = t55 * t28 - t41 * t29;
t31 = qJD(4) + (-pkin(3) * t47 + t50) * qJD(1);
t46 = cos(qJ(5));
t44 = sin(qJ(5));
t40 = qJD(3) + qJD(5);
t36 = t50 * qJD(1);
t33 = (t41 * t47 + t55 * t45) * qJD(1);
t32 = t41 * t54 - t55 * t53;
t24 = t32 * pkin(4) + t31;
t23 = -t44 * t32 + t46 * t33;
t22 = t46 * t32 + t44 * t33;
t19 = -t32 * pkin(7) + t21;
t18 = qJD(3) * pkin(4) - t33 * pkin(7) + t20;
t1 = [t57, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t42 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t48, t45 ^ 2 * t57, t45 * t48 * t47, t45 * t51, t47 * t51, qJD(3) ^ 2 / 0.2e1, -t36 * t53 + (-t45 * t35 + t39) * qJD(3), -t56 * qJD(3) + t36 * t54, t20 * qJD(3) + t31 * t32, -t21 * qJD(3) + t31 * t33, -t20 * t33 - t21 * t32, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t22, t23 * t40, -t22 * t40, t40 ^ 2 / 0.2e1, t24 * t22 + (t46 * t18 - t44 * t19) * t40, t24 * t23 - (t44 * t18 + t46 * t19) * t40;];
T_reg = t1;
