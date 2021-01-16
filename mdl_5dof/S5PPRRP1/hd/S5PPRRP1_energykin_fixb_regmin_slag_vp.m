% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:48:15
% EndTime: 2021-01-15 14:48:15
% DurationCPUTime: 0.05s
% Computational Cost: add. (61->23), mult. (161->59), div. (0->0), fcn. (98->6), ass. (0->23)
t18 = qJD(3) ^ 2;
t29 = t18 / 0.2e1;
t12 = sin(pkin(8));
t13 = cos(pkin(8));
t25 = qJD(1) * cos(qJ(3));
t26 = qJD(1) * sin(qJ(3));
t28 = t12 * t25 + t13 * t26;
t14 = sin(qJ(4));
t16 = cos(qJ(4));
t6 = qJD(3) * pkin(6) + t28;
t27 = t14 * qJD(2) + t16 * t6;
t24 = qJD(3) * t14;
t23 = qJD(3) * t16;
t22 = qJ(5) * qJD(3);
t21 = qJD(3) * qJD(4);
t20 = -t12 * t26 + t13 * t25;
t19 = qJD(1) ^ 2;
t11 = t16 * qJD(2);
t5 = -qJD(3) * pkin(3) - t20;
t3 = qJD(5) + (-pkin(4) * t16 - pkin(3)) * qJD(3) - t20;
t2 = t16 * t22 + t27;
t1 = qJD(4) * pkin(4) + t11 + (-t6 - t22) * t14;
t4 = [t19 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t12 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1) * t19, t29, t20 * qJD(3), -t28 * qJD(3), t14 ^ 2 * t29, t14 * t18 * t16, t14 * t21, t16 * t21, qJD(4) ^ 2 / 0.2e1, (-t14 * t6 + t11) * qJD(4) - t5 * t23, -t27 * qJD(4) + t5 * t24, t1 * qJD(4) - t3 * t23, -t2 * qJD(4) + t3 * t24, (-t1 * t14 + t16 * t2) * qJD(3), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t4;
