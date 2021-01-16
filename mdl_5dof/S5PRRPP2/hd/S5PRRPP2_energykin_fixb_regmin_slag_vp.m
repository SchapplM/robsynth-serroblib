% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPP2
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
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:31
% EndTime: 2021-01-15 15:32:31
% DurationCPUTime: 0.06s
% Computational Cost: add. (136->29), mult. (309->70), div. (0->0), fcn. (179->6), ass. (0->28)
t33 = qJD(2) ^ 2;
t41 = t33 / 0.2e1;
t29 = sin(qJ(3));
t30 = sin(qJ(2));
t24 = qJD(2) * pkin(6) + t30 * qJD(1);
t34 = qJ(4) * qJD(2) + t24;
t19 = qJD(3) * pkin(3) - t34 * t29;
t31 = cos(qJ(3));
t20 = t34 * t31;
t27 = sin(pkin(8));
t28 = cos(pkin(8));
t16 = t27 * t19 + t28 * t20;
t40 = qJD(2) * t29;
t39 = qJD(2) * t31;
t38 = qJD(3) * t24;
t32 = cos(qJ(2));
t37 = t32 * qJD(1);
t36 = qJD(1) * qJD(2);
t35 = qJD(2) * qJD(3);
t15 = t28 * t19 - t27 * t20;
t23 = -t37 + qJD(4) + (-pkin(3) * t31 - pkin(2)) * qJD(2);
t25 = -qJD(2) * pkin(2) - t37;
t22 = (t27 * t31 + t28 * t29) * qJD(2);
t21 = t27 * t40 - t28 * t39;
t14 = qJD(3) * qJ(5) + t16;
t13 = t21 * pkin(4) - t22 * qJ(5) + t23;
t12 = -qJD(3) * pkin(4) + qJD(5) - t15;
t1 = [qJD(1) ^ 2 / 0.2e1, t41, t32 * t36, -t30 * t36, t29 ^ 2 * t41, t29 * t33 * t31, t29 * t35, t31 * t35, qJD(3) ^ 2 / 0.2e1, -t25 * t39 - t29 * t38, t25 * t40 - t31 * t38, t15 * qJD(3) + t23 * t21, -t16 * qJD(3) + t23 * t22, -t15 * t22 - t16 * t21, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1, -t12 * qJD(3) + t13 * t21, t12 * t22 - t14 * t21, t14 * qJD(3) - t13 * t22, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1;];
T_reg = t1;
