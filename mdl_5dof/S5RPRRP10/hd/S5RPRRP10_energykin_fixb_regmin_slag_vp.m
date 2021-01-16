% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:14:47
% EndTime: 2021-01-15 19:14:47
% DurationCPUTime: 0.08s
% Computational Cost: add. (215->39), mult. (551->80), div. (0->0), fcn. (373->6), ass. (0->32)
t51 = cos(qJ(4));
t50 = pkin(6) + qJ(2);
t40 = sin(qJ(3));
t41 = cos(qJ(3));
t38 = cos(pkin(8));
t46 = qJD(1) * t38;
t37 = sin(pkin(8));
t47 = qJD(1) * t37;
t27 = t40 * t47 - t41 * t46;
t28 = (t37 * t41 + t38 * t40) * qJD(1);
t31 = qJD(2) + (-pkin(2) * t38 - pkin(1)) * qJD(1);
t17 = t27 * pkin(3) - t28 * pkin(7) + t31;
t29 = t50 * t47;
t30 = t50 * t46;
t48 = -t40 * t29 + t41 * t30;
t20 = qJD(3) * pkin(7) + t48;
t39 = sin(qJ(4));
t49 = t39 * t17 + t51 * t20;
t45 = t51 * t17 - t39 * t20;
t44 = -t41 * t29 - t40 * t30;
t19 = -qJD(3) * pkin(3) - t44;
t42 = qJD(1) ^ 2;
t36 = t38 ^ 2;
t35 = t37 ^ 2;
t33 = -qJD(1) * pkin(1) + qJD(2);
t23 = qJD(4) + t27;
t22 = t39 * qJD(3) + t51 * t28;
t21 = -t51 * qJD(3) + t39 * t28;
t14 = t21 * pkin(4) + qJD(5) + t19;
t13 = -t21 * qJ(5) + t49;
t12 = t23 * pkin(4) - t22 * qJ(5) + t45;
t1 = [t42 / 0.2e1, 0, 0, -t33 * t46, (t35 + t36) * t42 * qJ(2), t33 ^ 2 / 0.2e1 + (t36 / 0.2e1 + t35 / 0.2e1) * qJ(2) ^ 2 * t42, t28 ^ 2 / 0.2e1, -t28 * t27, t28 * qJD(3), -t27 * qJD(3), qJD(3) ^ 2 / 0.2e1, t44 * qJD(3) + t31 * t27, -t48 * qJD(3) + t31 * t28, t22 ^ 2 / 0.2e1, -t22 * t21, t22 * t23, -t21 * t23, t23 ^ 2 / 0.2e1, t19 * t21 + t45 * t23, t19 * t22 - t49 * t23, t12 * t23 + t14 * t21, -t13 * t23 + t14 * t22, -t12 * t22 - t13 * t21, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1;];
T_reg = t1;
