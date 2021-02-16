% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:50
% EndTime: 2021-01-15 15:04:51
% DurationCPUTime: 0.12s
% Computational Cost: add. (89->28), mult. (212->65), div. (0->0), fcn. (118->4), ass. (0->25)
t22 = sin(pkin(8));
t26 = qJD(2) ^ 2;
t36 = t22 ^ 2 * t26;
t23 = cos(pkin(8));
t12 = qJD(3) + (-pkin(3) * t23 - pkin(6) * t22 - pkin(2)) * qJD(2);
t32 = qJ(3) * qJD(2);
t16 = qJD(1) * t22 + t23 * t32;
t24 = sin(qJ(4));
t25 = cos(qJ(4));
t35 = t24 * t12 + t25 * t16;
t34 = qJD(2) * t22;
t33 = qJD(2) * t23;
t31 = t24 * t34;
t30 = t25 * t34;
t20 = t23 * qJD(1);
t14 = t22 * t32 - t20;
t29 = t14 * t34;
t28 = qJ(5) * t34;
t27 = t25 * t12 - t16 * t24;
t19 = -qJD(2) * pkin(2) + qJD(3);
t17 = -qJD(4) + t33;
t9 = qJD(5) - t20 + (pkin(4) * t24 + qJ(3)) * t34;
t8 = -t24 * t28 + t35;
t7 = -pkin(4) * t17 - t25 * t28 + t27;
t1 = [qJD(1) ^ 2 / 0.2e1, t26 / 0.2e1, 0, 0, -t19 * t33, (t14 * t22 + t16 * t23) * qJD(2), t16 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t25 ^ 2 * t36 / 0.2e1, -t25 * t24 * t36, -t17 * t30, t17 * t31, t17 ^ 2 / 0.2e1, -t27 * t17 + t24 * t29, t35 * t17 + t25 * t29, -t17 * t7 + t9 * t31, t17 * t8 + t9 * t30, (-t24 * t8 - t25 * t7) * t34, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg = t1;
