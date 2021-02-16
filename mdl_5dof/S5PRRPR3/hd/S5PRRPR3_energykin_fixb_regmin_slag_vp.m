% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:08
% EndTime: 2021-01-15 15:42:08
% DurationCPUTime: 0.08s
% Computational Cost: add. (127->32), mult. (321->73), div. (0->0), fcn. (210->6), ass. (0->28)
t39 = qJD(2) ^ 2;
t46 = t39 / 0.2e1;
t38 = cos(qJ(3));
t45 = t38 * t39;
t32 = t38 * qJD(1);
t36 = sin(qJ(3));
t42 = qJD(2) * t36;
t23 = qJD(3) * pkin(3) + t32 + (-pkin(6) - qJ(4)) * t42;
t41 = qJD(2) * t38;
t44 = pkin(6) * t41 + t36 * qJD(1);
t25 = qJ(4) * t41 + t44;
t34 = sin(pkin(9));
t43 = cos(pkin(9));
t16 = t34 * t23 + t43 * t25;
t40 = qJD(2) * qJD(3);
t15 = t43 * t23 - t34 * t25;
t28 = qJD(4) + (-pkin(3) * t38 - pkin(2)) * qJD(2);
t37 = cos(qJ(5));
t35 = sin(qJ(5));
t33 = qJD(3) + qJD(5);
t27 = (t34 * t38 + t43 * t36) * qJD(2);
t26 = t34 * t42 - t43 * t41;
t19 = t26 * pkin(4) + t28;
t18 = -t35 * t26 + t37 * t27;
t17 = t37 * t26 + t35 * t27;
t14 = -t26 * pkin(7) + t16;
t13 = qJD(3) * pkin(4) - t27 * pkin(7) + t15;
t1 = [qJD(1) ^ 2 / 0.2e1, t46, 0, 0, t36 ^ 2 * t46, t36 * t45, t36 * t40, t38 * t40, qJD(3) ^ 2 / 0.2e1, pkin(2) * t45 + (-pkin(6) * t42 + t32) * qJD(3), -t39 * pkin(2) * t36 - t44 * qJD(3), t15 * qJD(3) + t28 * t26, -t16 * qJD(3) + t28 * t27, -t15 * t27 - t16 * t26, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t17, t18 * t33, -t17 * t33, t33 ^ 2 / 0.2e1, t19 * t17 + (t37 * t13 - t35 * t14) * t33, t19 * t18 - (t35 * t13 + t37 * t14) * t33;];
T_reg = t1;
