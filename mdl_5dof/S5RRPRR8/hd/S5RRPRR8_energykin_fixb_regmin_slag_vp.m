% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:10
% EndTime: 2021-01-15 21:35:10
% DurationCPUTime: 0.10s
% Computational Cost: add. (252->41), mult. (650->90), div. (0->0), fcn. (464->8), ass. (0->39)
t54 = qJD(1) ^ 2;
t65 = t54 / 0.2e1;
t53 = cos(qJ(2));
t64 = t53 * t54;
t63 = pkin(6) + qJ(3);
t50 = sin(qJ(2));
t60 = qJD(1) * t50;
t41 = qJD(2) * pkin(2) - t63 * t60;
t59 = qJD(1) * t53;
t42 = t63 * t59;
t47 = sin(pkin(9));
t61 = cos(pkin(9));
t32 = t61 * t41 - t47 * t42;
t39 = (t47 * t53 + t61 * t50) * qJD(1);
t25 = qJD(2) * pkin(3) - t39 * pkin(7) + t32;
t33 = t47 * t41 + t61 * t42;
t38 = t47 * t60 - t61 * t59;
t26 = -t38 * pkin(7) + t33;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t62 = t49 * t25 + t52 * t26;
t58 = qJD(1) * qJD(2);
t57 = t50 * t58;
t56 = t53 * t58;
t30 = t52 * t38 + t49 * t39;
t55 = t52 * t25 - t49 * t26;
t43 = qJD(3) + (-pkin(2) * t53 - pkin(1)) * qJD(1);
t34 = t38 * pkin(3) + t43;
t51 = cos(qJ(5));
t48 = sin(qJ(5));
t46 = qJD(2) + qJD(4);
t31 = -t49 * t38 + t52 * t39;
t29 = qJD(5) + t30;
t28 = t51 * t31 + t48 * t46;
t27 = t48 * t31 - t51 * t46;
t22 = t30 * pkin(4) - t31 * pkin(8) + t34;
t21 = t46 * pkin(8) + t62;
t20 = -t46 * pkin(4) - t55;
t1 = [t65, 0, 0, t50 ^ 2 * t65, t50 * t64, t57, t56, qJD(2) ^ 2 / 0.2e1, pkin(1) * t64 - pkin(6) * t57, -t54 * pkin(1) * t50 - pkin(6) * t56, t32 * qJD(2) + t43 * t38, -t33 * qJD(2) + t43 * t39, -t32 * t39 - t33 * t38, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1, t31 ^ 2 / 0.2e1, -t31 * t30, t31 * t46, -t30 * t46, t46 ^ 2 / 0.2e1, t34 * t30 + t55 * t46, t34 * t31 - t62 * t46, t28 ^ 2 / 0.2e1, -t28 * t27, t28 * t29, -t27 * t29, t29 ^ 2 / 0.2e1, (-t48 * t21 + t51 * t22) * t29 + t20 * t27, -(t51 * t21 + t48 * t22) * t29 + t20 * t28;];
T_reg = t1;
