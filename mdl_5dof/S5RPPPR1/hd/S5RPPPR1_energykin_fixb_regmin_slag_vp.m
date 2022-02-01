% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:43
% EndTime: 2022-01-20 09:12:43
% DurationCPUTime: 0.19s
% Computational Cost: add. (117->34), mult. (291->79), div. (0->0), fcn. (169->8), ass. (0->29)
t37 = sin(pkin(8));
t39 = cos(pkin(9));
t50 = t37 * t39;
t40 = cos(pkin(8));
t41 = cos(pkin(7));
t46 = -pkin(1) * t41 - pkin(2);
t24 = qJD(3) + (-pkin(3) * t40 - qJ(4) * t37 + t46) * qJD(1);
t38 = sin(pkin(7));
t33 = (pkin(1) * t38 + qJ(3)) * qJD(1);
t30 = t37 * qJD(2) + t40 * t33;
t36 = sin(pkin(9));
t20 = t36 * t24 + t39 * t30;
t49 = qJD(1) * t37;
t48 = t40 * qJD(1);
t47 = t36 * t49;
t19 = t39 * t24 - t36 * t30;
t29 = t40 * qJD(2) - t37 * t33;
t28 = qJD(4) - t29;
t44 = qJD(1) ^ 2;
t43 = cos(qJ(5));
t42 = sin(qJ(5));
t34 = -qJD(5) + t48;
t32 = t46 * qJD(1) + qJD(3);
t27 = (-t36 * t42 + t39 * t43) * t49;
t26 = (t36 * t43 + t39 * t42) * t49;
t21 = pkin(4) * t47 + t28;
t18 = -pkin(6) * t47 + t20;
t17 = (-pkin(4) * t40 - pkin(6) * t50) * qJD(1) + t19;
t1 = [t44 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t38 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t44, -t32 * t48, (-t29 * t37 + t30 * t40) * qJD(1), t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, (t28 * t36 * t37 - t19 * t40) * qJD(1), (t20 * t40 + t28 * t50) * qJD(1), t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t26, -t27 * t34, t26 * t34, t34 ^ 2 / 0.2e1, -(t43 * t17 - t42 * t18) * t34 + t21 * t26, (t42 * t17 + t43 * t18) * t34 + t21 * t27;];
T_reg = t1;
