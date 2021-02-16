% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:22:38
% EndTime: 2021-01-15 15:22:38
% DurationCPUTime: 0.06s
% Computational Cost: add. (127->29), mult. (297->66), div. (0->0), fcn. (169->4), ass. (0->23)
t32 = qJD(2) ^ 2;
t38 = t32 / 0.2e1;
t31 = cos(qJ(3));
t37 = t31 * t32;
t27 = t31 * qJD(1);
t30 = sin(qJ(3));
t35 = qJD(2) * t30;
t19 = qJD(3) * pkin(3) + t27 + (-pkin(6) - qJ(4)) * t35;
t34 = qJD(2) * t31;
t36 = pkin(6) * t34 + t30 * qJD(1);
t20 = qJ(4) * t34 + t36;
t28 = sin(pkin(8));
t29 = cos(pkin(8));
t16 = t28 * t19 + t29 * t20;
t33 = qJD(2) * qJD(3);
t15 = t29 * t19 - t28 * t20;
t23 = qJD(4) + (-pkin(3) * t31 - pkin(2)) * qJD(2);
t22 = (t28 * t31 + t29 * t30) * qJD(2);
t21 = t28 * t35 - t29 * t34;
t14 = t21 * pkin(4) - t22 * qJ(5) + t23;
t13 = qJD(3) * qJ(5) + t16;
t12 = -qJD(3) * pkin(4) + qJD(5) - t15;
t1 = [qJD(1) ^ 2 / 0.2e1, t38, 0, 0, t30 ^ 2 * t38, t30 * t37, t30 * t33, t31 * t33, qJD(3) ^ 2 / 0.2e1, pkin(2) * t37 + (-pkin(6) * t35 + t27) * qJD(3), -t32 * pkin(2) * t30 - t36 * qJD(3), t15 * qJD(3) + t23 * t21, -t16 * qJD(3) + t23 * t22, -t15 * t22 - t16 * t21, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, -t12 * qJD(3) + t14 * t21, t12 * t22 - t13 * t21, t13 * qJD(3) - t14 * t22, t13 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1;];
T_reg = t1;
