% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:39
% EndTime: 2021-01-15 16:14:39
% DurationCPUTime: 0.05s
% Computational Cost: add. (79->20), mult. (123->51), div. (0->0), fcn. (51->4), ass. (0->20)
t12 = qJD(2) + qJD(3);
t11 = t12 ^ 2;
t25 = t11 / 0.2e1;
t13 = sin(qJ(4));
t15 = cos(qJ(4));
t21 = pkin(2) * qJD(2);
t18 = sin(qJ(3)) * t21;
t7 = t12 * pkin(7) + t18;
t24 = t13 * qJD(1) + t15 * t7;
t23 = t12 * t13;
t22 = t12 * t15;
t20 = qJ(5) * t12;
t19 = qJD(4) * t12;
t17 = cos(qJ(3)) * t21;
t10 = t15 * qJD(1);
t8 = -t12 * pkin(3) - t17;
t5 = -t17 + qJD(5) + (-pkin(4) * t15 - pkin(3)) * t12;
t4 = t15 * t20 + t24;
t3 = qJD(4) * pkin(4) + t10 + (-t7 - t20) * t13;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t25, t12 * t17, -t12 * t18, t13 ^ 2 * t25, t13 * t11 * t15, t13 * t19, t15 * t19, qJD(4) ^ 2 / 0.2e1, -t8 * t22 + (-t13 * t7 + t10) * qJD(4), -t24 * qJD(4) + t8 * t23, t3 * qJD(4) - t5 * t22, -t4 * qJD(4) + t5 * t23, (-t13 * t3 + t15 * t4) * t12, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t1;
