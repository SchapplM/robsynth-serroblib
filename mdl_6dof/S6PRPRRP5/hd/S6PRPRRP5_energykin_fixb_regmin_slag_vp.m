% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:51:49
% EndTime: 2021-01-16 01:51:49
% DurationCPUTime: 0.08s
% Computational Cost: add. (173->38), mult. (352->80), div. (0->0), fcn. (217->8), ass. (0->36)
t36 = qJD(2) ^ 2;
t53 = t36 / 0.2e1;
t52 = qJD(1) ^ 2 / 0.2e1;
t51 = cos(qJ(5));
t35 = cos(qJ(2));
t48 = qJD(1) * sin(pkin(6));
t38 = -t35 * t48 + qJD(3);
t19 = (-pkin(8) - pkin(2)) * qJD(2) + t38;
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t30 = cos(pkin(6));
t47 = qJD(1) * t30;
t49 = t32 * t19 + t34 * t47;
t14 = qJD(4) * pkin(9) + t49;
t33 = sin(qJ(2));
t42 = t33 * t48;
t17 = t42 + (pkin(4) * t32 - pkin(9) * t34 + qJ(3)) * qJD(2);
t31 = sin(qJ(5));
t50 = t51 * t14 + t31 * t17;
t46 = qJD(2) * t34;
t23 = qJD(2) * qJ(3) + t42;
t45 = t23 * qJD(2);
t44 = t32 * qJD(2);
t43 = qJD(2) * qJD(4);
t41 = qJD(2) * t48;
t40 = -t31 * t14 + t51 * t17;
t39 = t34 * t19 - t32 * t47;
t13 = -qJD(4) * pkin(4) - t39;
t27 = qJD(5) + t44;
t22 = t31 * qJD(4) + t51 * t46;
t21 = -t51 * qJD(4) + t31 * t46;
t20 = -qJD(2) * pkin(2) + t38;
t11 = t21 * pkin(5) + qJD(6) + t13;
t10 = -t21 * qJ(6) + t50;
t9 = t27 * pkin(5) - t22 * qJ(6) + t40;
t1 = [t52, t53, t35 * t41, -t33 * t41, t20 * qJD(2), t45, t30 ^ 2 * t52 + t23 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t34 ^ 2 * t53, -t34 * t36 * t32, t34 * t43, -t32 * t43, qJD(4) ^ 2 / 0.2e1, t39 * qJD(4) + t23 * t44, -t49 * qJD(4) + t34 * t45, t22 ^ 2 / 0.2e1, -t22 * t21, t22 * t27, -t21 * t27, t27 ^ 2 / 0.2e1, t13 * t21 + t40 * t27, t13 * t22 - t50 * t27, t11 * t21 + t9 * t27, -t10 * t27 + t11 * t22, -t10 * t21 - t9 * t22, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1;];
T_reg = t1;
