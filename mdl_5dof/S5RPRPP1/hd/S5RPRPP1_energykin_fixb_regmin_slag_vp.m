% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:02
% EndTime: 2021-01-15 11:14:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (163->33), mult. (377->77), div. (0->0), fcn. (205->6), ass. (0->28)
t41 = qJD(1) ^ 2;
t49 = t41 / 0.2e1;
t36 = sin(pkin(7));
t30 = (pkin(1) * t36 + pkin(6)) * qJD(1);
t40 = cos(qJ(3));
t34 = t40 * qJD(2);
t39 = sin(qJ(3));
t45 = qJ(4) * qJD(1);
t24 = qJD(3) * pkin(3) + t34 + (-t30 - t45) * t39;
t48 = t39 * qJD(2) + t40 * t30;
t25 = t40 * t45 + t48;
t35 = sin(pkin(8));
t37 = cos(pkin(8));
t20 = t35 * t24 + t37 * t25;
t47 = qJD(1) * t39;
t46 = qJD(1) * t40;
t44 = qJD(1) * qJD(3);
t38 = cos(pkin(7));
t43 = -pkin(1) * t38 - pkin(2);
t19 = t37 * t24 - t35 * t25;
t26 = qJD(4) + (-pkin(3) * t40 + t43) * qJD(1);
t31 = t43 * qJD(1);
t28 = (t35 * t40 + t37 * t39) * qJD(1);
t27 = t35 * t47 - t37 * t46;
t21 = t27 * pkin(4) - t28 * qJ(5) + t26;
t18 = qJD(3) * qJ(5) + t20;
t17 = -qJD(3) * pkin(4) + qJD(5) - t19;
t1 = [t49, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t36 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t41, t39 ^ 2 * t49, t39 * t41 * t40, t39 * t44, t40 * t44, qJD(3) ^ 2 / 0.2e1, -t31 * t46 + (-t39 * t30 + t34) * qJD(3), -t48 * qJD(3) + t31 * t47, qJD(3) * t19 + t26 * t27, -qJD(3) * t20 + t26 * t28, -t19 * t28 - t20 * t27, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, -qJD(3) * t17 + t21 * t27, t17 * t28 - t18 * t27, qJD(3) * t18 - t21 * t28, t18 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1;];
T_reg = t1;
