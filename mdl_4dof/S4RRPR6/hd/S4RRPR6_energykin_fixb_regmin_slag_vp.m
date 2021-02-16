% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:20
% EndTime: 2021-01-15 10:46:20
% DurationCPUTime: 0.07s
% Computational Cost: add. (105->28), mult. (297->67), div. (0->0), fcn. (188->6), ass. (0->29)
t36 = qJD(1) ^ 2;
t45 = t36 / 0.2e1;
t35 = cos(qJ(2));
t44 = t35 * t36;
t43 = pkin(5) + qJ(3);
t33 = sin(qJ(2));
t41 = qJD(1) * t33;
t26 = qJD(2) * pkin(2) - t43 * t41;
t40 = qJD(1) * t35;
t27 = t43 * t40;
t31 = sin(pkin(7));
t42 = cos(pkin(7));
t18 = t31 * t26 + t42 * t27;
t39 = qJD(1) * qJD(2);
t38 = t33 * t39;
t37 = t35 * t39;
t17 = t42 * t26 - t31 * t27;
t28 = qJD(3) + (-pkin(2) * t35 - pkin(1)) * qJD(1);
t34 = cos(qJ(4));
t32 = sin(qJ(4));
t30 = qJD(2) + qJD(4);
t24 = (t31 * t35 + t42 * t33) * qJD(1);
t23 = t31 * t41 - t42 * t40;
t19 = t23 * pkin(3) + t28;
t16 = -t32 * t23 + t34 * t24;
t15 = t34 * t23 + t32 * t24;
t14 = -t23 * pkin(6) + t18;
t13 = qJD(2) * pkin(3) - t24 * pkin(6) + t17;
t1 = [t45, 0, 0, t33 ^ 2 * t45, t33 * t44, t38, t37, qJD(2) ^ 2 / 0.2e1, pkin(1) * t44 - pkin(5) * t38, -t36 * pkin(1) * t33 - pkin(5) * t37, t17 * qJD(2) + t28 * t23, -t18 * qJD(2) + t28 * t24, -t17 * t24 - t18 * t23, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t16 ^ 2 / 0.2e1, -t16 * t15, t16 * t30, -t15 * t30, t30 ^ 2 / 0.2e1, t19 * t15 + (t34 * t13 - t32 * t14) * t30, t19 * t16 - (t32 * t13 + t34 * t14) * t30;];
T_reg = t1;
