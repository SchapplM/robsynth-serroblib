% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:54
% EndTime: 2022-01-23 09:20:54
% DurationCPUTime: 0.06s
% Computational Cost: add. (107->26), mult. (202->63), div. (0->0), fcn. (103->8), ass. (0->28)
t25 = qJD(1) + qJD(3);
t23 = t25 ^ 2;
t26 = sin(pkin(9));
t44 = t23 * t26 ^ 2;
t43 = t25 * t26;
t28 = cos(pkin(9));
t42 = t28 * t25;
t29 = cos(pkin(8));
t17 = (pkin(1) * t29 + pkin(2)) * qJD(1);
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t27 = sin(pkin(8));
t38 = pkin(1) * qJD(1) * t27;
t41 = t31 * t17 + t33 * t38;
t30 = sin(qJ(5));
t40 = t30 * t43;
t32 = cos(qJ(5));
t39 = t32 * t43;
t37 = t33 * t17 - t31 * t38;
t36 = qJD(4) - t37;
t34 = qJD(1) ^ 2;
t18 = -qJD(5) + t42;
t15 = t25 * qJ(4) + t41;
t14 = -t25 * pkin(3) + t36;
t13 = t26 * qJD(2) + t28 * t15;
t11 = -t28 * qJD(2) + t26 * t15;
t10 = (-pkin(4) * t28 - pkin(7) * t26 - pkin(3)) * t25 + t36;
t1 = [t34 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t27 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t34, t23 / 0.2e1, t37 * t25, -t41 * t25, -t14 * t42, (t11 * t26 + t13 * t28) * t25, t13 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t32 ^ 2 * t44 / 0.2e1, -t32 * t30 * t44, -t18 * t39, t18 * t40, t18 ^ 2 / 0.2e1, -(t32 * t10 - t30 * t13) * t18 + t11 * t40, (t30 * t10 + t32 * t13) * t18 + t11 * t39;];
T_reg = t1;
