% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:27
% EndTime: 2021-01-15 14:56:27
% DurationCPUTime: 0.04s
% Computational Cost: add. (45->20), mult. (110->50), div. (0->0), fcn. (51->4), ass. (0->21)
t12 = qJD(3) ^ 2;
t21 = t12 / 0.2e1;
t8 = sin(qJ(4));
t20 = qJD(3) * t8;
t10 = cos(qJ(4));
t19 = qJD(3) * t10;
t18 = t10 * qJD(1);
t11 = cos(qJ(3));
t17 = t11 * qJD(2);
t16 = qJ(5) * qJD(3);
t15 = qJD(2) * qJD(3);
t14 = qJD(3) * qJD(4);
t9 = sin(qJ(3));
t5 = qJD(3) * pkin(6) + t9 * qJD(2);
t13 = -t8 * qJD(1) + t10 * t5;
t7 = qJD(1) ^ 2 / 0.2e1;
t6 = -qJD(3) * pkin(3) - t17;
t3 = -t17 + qJD(5) + (-pkin(4) * t10 - pkin(3)) * qJD(3);
t2 = t10 * t16 + t13;
t1 = qJD(4) * pkin(4) - t18 + (-t5 - t16) * t8;
t4 = [t7, t7 + qJD(2) ^ 2 / 0.2e1, t21, t11 * t15, -t9 * t15, t8 ^ 2 * t21, t8 * t12 * t10, t8 * t14, t10 * t14, qJD(4) ^ 2 / 0.2e1, -t6 * t19 + (-t8 * t5 - t18) * qJD(4), -t13 * qJD(4) + t6 * t20, t1 * qJD(4) - t3 * t19, -t2 * qJD(4) + t3 * t20, (-t1 * t8 + t10 * t2) * qJD(3), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t4;
