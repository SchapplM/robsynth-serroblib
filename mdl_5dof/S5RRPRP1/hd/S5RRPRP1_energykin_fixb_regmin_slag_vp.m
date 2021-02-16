% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:08:51
% EndTime: 2021-01-15 20:08:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (126->25), mult. (206->60), div. (0->0), fcn. (98->6), ass. (0->25)
t22 = qJD(1) + qJD(2);
t21 = t22 ^ 2;
t37 = t21 / 0.2e1;
t25 = sin(qJ(4));
t36 = t22 * t25;
t27 = cos(qJ(4));
t35 = t22 * t27;
t33 = pkin(1) * qJD(1);
t29 = cos(qJ(2)) * t33;
t16 = t22 * pkin(2) + t29;
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t30 = sin(qJ(2)) * t33;
t14 = t23 * t16 + t24 * t30;
t12 = t22 * pkin(7) + t14;
t34 = t25 * qJD(3) + t27 * t12;
t32 = qJ(5) * t22;
t31 = qJD(4) * t22;
t13 = t24 * t16 - t23 * t30;
t20 = t27 * qJD(3);
t11 = -t22 * pkin(3) - t13;
t9 = qJD(5) + (-pkin(4) * t27 - pkin(3)) * t22 - t13;
t8 = t27 * t32 + t34;
t7 = qJD(4) * pkin(4) + t20 + (-t12 - t32) * t25;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t37, t22 * t29, -t22 * t30, t14 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t25 ^ 2 * t37, t25 * t21 * t27, t25 * t31, t27 * t31, qJD(4) ^ 2 / 0.2e1, -t11 * t35 + (-t25 * t12 + t20) * qJD(4), -t34 * qJD(4) + t11 * t36, t7 * qJD(4) - t9 * t35, -t8 * qJD(4) + t9 * t36, (-t25 * t7 + t27 * t8) * t22, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg = t1;
