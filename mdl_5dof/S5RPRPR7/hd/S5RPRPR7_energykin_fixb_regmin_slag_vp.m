% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR7
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
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:06:25
% EndTime: 2021-01-15 12:06:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (154->36), mult. (375->84), div. (0->0), fcn. (225->8), ass. (0->33)
t57 = qJD(1) ^ 2;
t65 = t57 / 0.2e1;
t50 = sin(pkin(8));
t43 = (pkin(1) * t50 + pkin(6)) * qJD(1);
t56 = cos(qJ(3));
t48 = t56 * qJD(2);
t54 = sin(qJ(3));
t61 = qJ(4) * qJD(1);
t34 = qJD(3) * pkin(3) + t48 + (-t43 - t61) * t54;
t64 = t54 * qJD(2) + t56 * t43;
t37 = t56 * t61 + t64;
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t30 = t49 * t34 + t51 * t37;
t63 = qJD(1) * t54;
t62 = qJD(1) * t56;
t60 = qJD(1) * qJD(3);
t52 = cos(pkin(8));
t59 = -pkin(1) * t52 - pkin(2);
t40 = t49 * t63 - t51 * t62;
t29 = t51 * t34 - t49 * t37;
t39 = qJD(4) + (-pkin(3) * t56 + t59) * qJD(1);
t55 = cos(qJ(5));
t53 = sin(qJ(5));
t44 = t59 * qJD(1);
t41 = (t49 * t56 + t51 * t54) * qJD(1);
t38 = qJD(5) + t40;
t36 = t53 * qJD(3) + t55 * t41;
t35 = -t55 * qJD(3) + t53 * t41;
t31 = t40 * pkin(4) - t41 * pkin(7) + t39;
t28 = qJD(3) * pkin(7) + t30;
t27 = -qJD(3) * pkin(4) - t29;
t1 = [t65, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t50 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t57, t54 ^ 2 * t65, t54 * t57 * t56, t54 * t60, t56 * t60, qJD(3) ^ 2 / 0.2e1, -t44 * t62 + (-t54 * t43 + t48) * qJD(3), -t64 * qJD(3) + t44 * t63, t29 * qJD(3) + t39 * t40, -t30 * qJD(3) + t39 * t41, -t29 * t41 - t30 * t40, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1, t36 ^ 2 / 0.2e1, -t36 * t35, t36 * t38, -t35 * t38, t38 ^ 2 / 0.2e1, (-t53 * t28 + t55 * t31) * t38 + t27 * t35, -(t55 * t28 + t53 * t31) * t38 + t27 * t36;];
T_reg = t1;
