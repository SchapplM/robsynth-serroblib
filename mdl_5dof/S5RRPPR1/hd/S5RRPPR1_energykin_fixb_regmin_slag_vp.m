% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:41
% EndTime: 2022-01-20 09:51:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (128->29), mult. (209->63), div. (0->0), fcn. (115->8), ass. (0->26)
t26 = qJD(1) + qJD(2);
t29 = cos(pkin(9));
t39 = t26 * t29;
t38 = pkin(1) * qJD(1);
t36 = cos(qJ(2)) * t38;
t20 = t26 * pkin(2) + t36;
t28 = sin(pkin(8));
t30 = cos(pkin(8));
t37 = sin(qJ(2)) * t38;
t16 = t28 * t20 + t30 * t37;
t14 = t26 * qJ(4) + t16;
t27 = sin(pkin(9));
t10 = t27 * qJD(3) + t29 * t14;
t15 = t30 * t20 - t28 * t37;
t35 = qJD(4) - t15;
t33 = cos(qJ(5));
t31 = sin(qJ(5));
t25 = t29 * qJD(3);
t18 = (t27 * t33 + t29 * t31) * t26;
t17 = t31 * t27 * t26 - t33 * t39;
t13 = -t26 * pkin(3) + t35;
t11 = (-pkin(4) * t29 - pkin(3)) * t26 + t35;
t9 = -t27 * t14 + t25;
t8 = pkin(7) * t39 + t10;
t7 = t25 + (-pkin(7) * t26 - t14) * t27;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t26 ^ 2 / 0.2e1, t26 * t36, -t26 * t37, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, -t13 * t39, (t10 * t29 - t27 * t9) * t26, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t17, t18 * qJD(5), -t17 * qJD(5), qJD(5) ^ 2 / 0.2e1, t11 * t17 + (-t31 * t8 + t33 * t7) * qJD(5), t11 * t18 - (t31 * t7 + t33 * t8) * qJD(5);];
T_reg = t1;
