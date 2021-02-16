% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP4
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
% Datum: 2021-01-15 12:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:56:34
% EndTime: 2021-01-15 12:56:34
% DurationCPUTime: 0.09s
% Computational Cost: add. (194->39), mult. (496->89), div. (0->0), fcn. (314->6), ass. (0->34)
t39 = sin(pkin(8));
t37 = t39 ^ 2;
t45 = qJD(1) ^ 2;
t57 = t37 * t45;
t40 = cos(pkin(8));
t27 = qJD(2) + (-pkin(2) * t40 - pkin(6) * t39 - pkin(1)) * qJD(1);
t44 = cos(qJ(3));
t26 = t44 * t27;
t52 = t40 * qJD(1);
t33 = -qJD(3) + t52;
t42 = sin(qJ(3));
t19 = -t33 * pkin(3) + t26 + (-pkin(7) * t39 * t44 - qJ(2) * t40 * t42) * qJD(1);
t53 = qJD(1) * t39;
t49 = t42 * t53;
t51 = qJ(2) * qJD(1);
t48 = t40 * t51;
t55 = t42 * t27 + t44 * t48;
t22 = -pkin(7) * t49 + t55;
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t56 = t41 * t19 + t43 * t22;
t28 = pkin(3) * t49 + t39 * t51;
t54 = qJ(2) * t45;
t50 = t37 * t54;
t47 = t43 * t19 - t41 * t22;
t38 = t40 ^ 2;
t36 = -qJD(1) * pkin(1) + qJD(2);
t30 = -qJD(4) + t33;
t24 = (-t41 * t42 + t43 * t44) * t53;
t23 = (t41 * t44 + t42 * t43) * t53;
t21 = t23 * pkin(4) + qJD(5) + t28;
t16 = -t23 * qJ(5) + t56;
t15 = -t30 * pkin(4) - t24 * qJ(5) + t47;
t1 = [t45 / 0.2e1, 0, 0, -t36 * t52, (t37 + t38) * t54, t36 ^ 2 / 0.2e1 + (t38 / 0.2e1 + t37 / 0.2e1) * qJ(2) ^ 2 * t45, t44 ^ 2 * t57 / 0.2e1, -t44 * t42 * t57, -t44 * t33 * t53, t33 * t49, t33 ^ 2 / 0.2e1, t42 * t50 - (-t42 * t48 + t26) * t33, t55 * t33 + t44 * t50, t24 ^ 2 / 0.2e1, -t24 * t23, -t24 * t30, t23 * t30, t30 ^ 2 / 0.2e1, t28 * t23 - t47 * t30, t28 * t24 + t56 * t30, -t15 * t30 + t21 * t23, t16 * t30 + t21 * t24, -t15 * t24 - t16 * t23, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1;];
T_reg = t1;
