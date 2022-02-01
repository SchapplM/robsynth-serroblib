% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:44
% EndTime: 2022-01-23 09:12:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (118->32), mult. (278->75), div. (0->0), fcn. (147->6), ass. (0->29)
t26 = sin(pkin(8));
t32 = qJD(1) ^ 2;
t44 = t26 ^ 2 * t32;
t28 = cos(pkin(8));
t29 = cos(pkin(7));
t36 = -pkin(1) * t29 - pkin(2);
t15 = qJD(3) + (-pkin(3) * t28 - pkin(6) * t26 + t36) * qJD(1);
t27 = sin(pkin(7));
t21 = (pkin(1) * t27 + qJ(3)) * qJD(1);
t19 = t26 * qJD(2) + t28 * t21;
t30 = sin(qJ(4));
t31 = cos(qJ(4));
t43 = t30 * t15 + t31 * t19;
t42 = qJD(1) * t26;
t41 = qJD(1) * t30;
t40 = t28 * qJD(1);
t39 = t26 * t41;
t38 = t31 * t42;
t24 = t28 * qJD(2);
t17 = t26 * t21 - t24;
t37 = t17 * t42;
t35 = qJ(5) * t42;
t34 = t31 * t15 - t30 * t19;
t22 = -qJD(4) + t40;
t20 = t36 * qJD(1) + qJD(3);
t12 = qJD(5) - t24 + (pkin(4) * t41 + t21) * t26;
t11 = -t30 * t35 + t43;
t10 = -t22 * pkin(4) - t31 * t35 + t34;
t1 = [t32 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t27 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t32, -t20 * t40, (t17 * t26 + t19 * t28) * qJD(1), t19 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t31 ^ 2 * t44 / 0.2e1, -t31 * t30 * t44, -t22 * t38, t22 * t39, t22 ^ 2 / 0.2e1, -t34 * t22 + t30 * t37, t43 * t22 + t31 * t37, -t10 * t22 + t12 * t39, t11 * t22 + t12 * t38, (-t10 * t31 - t11 * t30) * t42, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1;];
T_reg = t1;
