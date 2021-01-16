% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:56:37
% EndTime: 2021-01-15 10:56:37
% DurationCPUTime: 0.07s
% Computational Cost: add. (101->28), mult. (281->67), div. (0->0), fcn. (172->6), ass. (0->29)
t45 = qJD(1) ^ 2;
t53 = t45 / 0.2e1;
t44 = cos(qJ(2));
t52 = t44 * t45;
t51 = pkin(5) + qJ(3);
t42 = sin(qJ(2));
t50 = qJD(1) * t42;
t34 = qJD(2) * pkin(2) - t51 * t50;
t49 = qJD(1) * t44;
t35 = t51 * t49;
t39 = sin(pkin(7));
t40 = cos(pkin(7));
t26 = t39 * t34 + t40 * t35;
t48 = qJD(1) * qJD(2);
t47 = t42 * t48;
t46 = t44 * t48;
t31 = t39 * t50 - t40 * t49;
t25 = t40 * t34 - t39 * t35;
t36 = qJD(3) + (-pkin(2) * t44 - pkin(1)) * qJD(1);
t43 = cos(qJ(4));
t41 = sin(qJ(4));
t32 = (t39 * t44 + t40 * t42) * qJD(1);
t30 = qJD(4) + t31;
t28 = t41 * qJD(2) + t43 * t32;
t27 = -t43 * qJD(2) + t41 * t32;
t24 = qJD(2) * pkin(6) + t26;
t23 = -qJD(2) * pkin(3) - t25;
t22 = t31 * pkin(3) - t32 * pkin(6) + t36;
t1 = [t53, 0, 0, t42 ^ 2 * t53, t42 * t52, t47, t46, qJD(2) ^ 2 / 0.2e1, pkin(1) * t52 - pkin(5) * t47, -t45 * pkin(1) * t42 - pkin(5) * t46, t25 * qJD(2) + t36 * t31, -t26 * qJD(2) + t36 * t32, -t25 * t32 - t26 * t31, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t28 ^ 2 / 0.2e1, -t28 * t27, t28 * t30, -t27 * t30, t30 ^ 2 / 0.2e1, (t43 * t22 - t41 * t24) * t30 + t23 * t27, -(t41 * t22 + t43 * t24) * t30 + t23 * t28;];
T_reg = t1;
