% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:16
% EndTime: 2021-01-16 01:24:16
% DurationCPUTime: 0.09s
% Computational Cost: add. (204->37), mult. (461->83), div. (0->0), fcn. (321->10), ass. (0->36)
t50 = qJD(2) ^ 2;
t62 = t50 / 0.2e1;
t61 = cos(qJ(5));
t49 = cos(qJ(2));
t58 = qJD(1) * sin(pkin(6));
t34 = qJD(2) * pkin(2) + t49 * t58;
t42 = sin(pkin(11));
t44 = cos(pkin(11));
t47 = sin(qJ(2));
t54 = t47 * t58;
t30 = t42 * t34 + t44 * t54;
t28 = qJD(2) * pkin(8) + t30;
t38 = cos(pkin(6)) * qJD(1) + qJD(3);
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t59 = t48 * t28 + t46 * t38;
t21 = qJD(4) * pkin(9) + t59;
t29 = t44 * t34 - t42 * t54;
t24 = (-pkin(4) * t48 - pkin(9) * t46 - pkin(3)) * qJD(2) - t29;
t45 = sin(qJ(5));
t60 = t61 * t21 + t45 * t24;
t57 = qJD(2) * t46;
t56 = t48 * qJD(2);
t55 = qJD(2) * qJD(4);
t53 = qJD(2) * t58;
t52 = -t45 * t21 + t61 * t24;
t51 = -t46 * t28 + t48 * t38;
t20 = -qJD(4) * pkin(4) - t51;
t39 = -qJD(5) + t56;
t33 = t45 * qJD(4) + t61 * t57;
t32 = -t61 * qJD(4) + t45 * t57;
t27 = -qJD(2) * pkin(3) - t29;
t18 = t32 * pkin(5) + qJD(6) + t20;
t17 = -t32 * qJ(6) + t60;
t16 = -t39 * pkin(5) - t33 * qJ(6) + t52;
t1 = [qJD(1) ^ 2 / 0.2e1, t62, t49 * t53, -t47 * t53, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t46 ^ 2 * t62, t46 * t50 * t48, t46 * t55, t48 * t55, qJD(4) ^ 2 / 0.2e1, t51 * qJD(4) - t27 * t56, -t59 * qJD(4) + t27 * t57, t33 ^ 2 / 0.2e1, -t33 * t32, -t33 * t39, t32 * t39, t39 ^ 2 / 0.2e1, t20 * t32 - t52 * t39, t20 * t33 + t60 * t39, -t16 * t39 + t18 * t32, t17 * t39 + t18 * t33, -t16 * t33 - t17 * t32, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1;];
T_reg = t1;
