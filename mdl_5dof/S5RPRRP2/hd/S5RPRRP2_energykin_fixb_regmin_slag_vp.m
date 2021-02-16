% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:36:26
% EndTime: 2021-01-15 12:36:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (111->25), mult. (210->62), div. (0->0), fcn. (98->6), ass. (0->26)
t20 = qJD(1) + qJD(3);
t19 = t20 ^ 2;
t37 = t19 / 0.2e1;
t23 = sin(qJ(4));
t36 = t20 * t23;
t25 = cos(qJ(4));
t35 = t20 * t25;
t22 = cos(pkin(8));
t14 = (pkin(1) * t22 + pkin(2)) * qJD(1);
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t21 = sin(pkin(8));
t30 = pkin(1) * qJD(1) * t21;
t33 = t24 * t14 + t26 * t30;
t12 = t20 * pkin(7) + t33;
t34 = t23 * qJD(2) + t25 * t12;
t32 = qJ(5) * t20;
t31 = qJD(4) * t20;
t29 = t26 * t14 - t24 * t30;
t27 = qJD(1) ^ 2;
t18 = t25 * qJD(2);
t11 = -t20 * pkin(3) - t29;
t9 = qJD(5) + (-pkin(4) * t25 - pkin(3)) * t20 - t29;
t8 = t25 * t32 + t34;
t7 = qJD(4) * pkin(4) + t18 + (-t12 - t32) * t23;
t1 = [t27 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t21 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t27, t37, t29 * t20, -t33 * t20, t23 ^ 2 * t37, t23 * t19 * t25, t23 * t31, t25 * t31, qJD(4) ^ 2 / 0.2e1, -t11 * t35 + (-t23 * t12 + t18) * qJD(4), -t34 * qJD(4) + t11 * t36, t7 * qJD(4) - t9 * t35, -t8 * qJD(4) + t9 * t36, (-t23 * t7 + t25 * t8) * t20, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1;];
T_reg = t1;
