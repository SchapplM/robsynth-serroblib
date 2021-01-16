% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRP3
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
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:29
% EndTime: 2021-01-16 01:38:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (296->46), mult. (715->93), div. (0->0), fcn. (539->10), ass. (0->37)
t59 = cos(qJ(5));
t47 = sin(qJ(2));
t56 = qJD(1) * sin(pkin(6));
t36 = qJD(2) * qJ(3) + t47 * t56;
t43 = cos(pkin(11));
t55 = qJD(1) * cos(pkin(6));
t38 = t43 * t55;
t41 = sin(pkin(11));
t24 = t38 + (-pkin(8) * qJD(2) - t36) * t41;
t29 = t43 * t36 + t41 * t55;
t54 = qJD(2) * t43;
t25 = pkin(8) * t54 + t29;
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t57 = t46 * t24 + t48 * t25;
t17 = qJD(4) * pkin(9) + t57;
t49 = cos(qJ(2));
t50 = -t49 * t56 + qJD(3);
t30 = (-pkin(3) * t43 - pkin(2)) * qJD(2) + t50;
t32 = t46 * t41 * qJD(2) - t48 * t54;
t33 = (t41 * t48 + t43 * t46) * qJD(2);
t20 = t32 * pkin(4) - t33 * pkin(9) + t30;
t45 = sin(qJ(5));
t58 = t59 * t17 + t45 * t20;
t53 = qJD(2) * t56;
t52 = -t45 * t17 + t59 * t20;
t51 = t48 * t24 - t46 * t25;
t16 = -qJD(4) * pkin(4) - t51;
t35 = -qJD(2) * pkin(2) + t50;
t31 = qJD(5) + t32;
t28 = -t41 * t36 + t38;
t27 = t45 * qJD(4) + t59 * t33;
t26 = -t59 * qJD(4) + t45 * t33;
t14 = t26 * pkin(5) + qJD(6) + t16;
t13 = -t26 * qJ(6) + t58;
t12 = t31 * pkin(5) - t27 * qJ(6) + t52;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t49 * t53, -t47 * t53, -t35 * t54, (-t28 * t41 + t29 * t43) * qJD(2), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t32, t33 * qJD(4), -t32 * qJD(4), qJD(4) ^ 2 / 0.2e1, t51 * qJD(4) + t30 * t32, -t57 * qJD(4) + t30 * t33, t27 ^ 2 / 0.2e1, -t27 * t26, t27 * t31, -t26 * t31, t31 ^ 2 / 0.2e1, t16 * t26 + t52 * t31, t16 * t27 - t58 * t31, t12 * t31 + t14 * t26, -t13 * t31 + t14 * t27, -t12 * t27 - t13 * t26, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1;];
T_reg = t1;
