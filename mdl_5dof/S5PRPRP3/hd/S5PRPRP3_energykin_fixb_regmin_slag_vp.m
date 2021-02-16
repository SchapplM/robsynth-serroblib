% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRP3
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
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:13:59
% EndTime: 2021-01-15 15:13:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (76->24), mult. (174->60), div. (0->0), fcn. (98->6), ass. (0->25)
t27 = qJD(2) ^ 2;
t35 = t27 / 0.2e1;
t26 = cos(qJ(2));
t16 = qJD(2) * pkin(2) + t26 * qJD(1);
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t24 = sin(qJ(2));
t33 = qJD(1) * t24;
t14 = t21 * t16 + t22 * t33;
t12 = qJD(2) * pkin(6) + t14;
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t34 = t23 * qJD(3) + t25 * t12;
t32 = qJD(2) * t23;
t31 = qJD(2) * t25;
t30 = qJ(5) * qJD(2);
t29 = qJD(1) * qJD(2);
t28 = qJD(2) * qJD(4);
t13 = t22 * t16 - t21 * t33;
t20 = t25 * qJD(3);
t11 = -qJD(2) * pkin(3) - t13;
t9 = qJD(5) + (-pkin(4) * t25 - pkin(3)) * qJD(2) - t13;
t8 = t25 * t30 + t34;
t7 = qJD(4) * pkin(4) + t20 + (-t12 - t30) * t23;
t1 = [qJD(1) ^ 2 / 0.2e1, t35, t26 * t29, -t24 * t29, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t23 ^ 2 * t35, t23 * t27 * t25, t23 * t28, t25 * t28, qJD(4) ^ 2 / 0.2e1, (-t23 * t12 + t20) * qJD(4) - t11 * t31, -t34 * qJD(4) + t11 * t32, t7 * qJD(4) - t9 * t31, -t8 * qJD(4) + t9 * t32, (-t23 * t7 + t25 * t8) * qJD(2), t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg = t1;
