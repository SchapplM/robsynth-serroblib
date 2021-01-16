% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:00
% EndTime: 2021-01-15 11:24:00
% DurationCPUTime: 0.07s
% Computational Cost: add. (158->31), mult. (306->69), div. (0->0), fcn. (143->4), ass. (0->24)
t34 = qJD(1) ^ 2;
t40 = t34 / 0.2e1;
t33 = cos(qJ(3));
t25 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t35 = -qJ(4) * qJD(1) + t25;
t20 = qJD(3) * pkin(3) + t35 * t33;
t32 = sin(qJ(3));
t21 = t35 * t32;
t30 = sin(pkin(7));
t31 = cos(pkin(7));
t16 = t30 * t20 + t31 * t21;
t39 = t34 * qJ(2);
t38 = qJD(1) * t32;
t37 = qJD(3) * t25;
t36 = qJD(1) * qJD(3);
t24 = pkin(3) * t38 + qJD(1) * qJ(2) + qJD(4);
t15 = t31 * t20 - t30 * t21;
t28 = -qJD(1) * pkin(1) + qJD(2);
t23 = t31 * t33 * qJD(1) - t30 * t38;
t22 = (t30 * t33 + t31 * t32) * qJD(1);
t17 = t22 * pkin(4) - t23 * qJ(5) + t24;
t14 = qJD(3) * qJ(5) + t16;
t13 = -qJD(3) * pkin(4) + qJD(5) - t15;
t1 = [t40, 0, 0, t28 * qJD(1), t39, qJ(2) ^ 2 * t40 + t28 ^ 2 / 0.2e1, t33 ^ 2 * t40, -t33 * t34 * t32, t33 * t36, -t32 * t36, qJD(3) ^ 2 / 0.2e1, t32 * t39 + t33 * t37, -t32 * t37 + t33 * t39, t15 * qJD(3) + t24 * t22, -t16 * qJD(3) + t24 * t23, -t15 * t23 - t16 * t22, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, -t13 * qJD(3) + t17 * t22, t13 * t23 - t14 * t22, t14 * qJD(3) - t17 * t23, t14 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg = t1;
