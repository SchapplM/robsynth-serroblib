% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:17
% EndTime: 2021-01-15 21:13:17
% DurationCPUTime: 0.07s
% Computational Cost: add. (106->32), mult. (270->77), div. (0->0), fcn. (163->6), ass. (0->30)
t38 = qJD(1) ^ 2;
t45 = t38 / 0.2e1;
t37 = cos(qJ(2));
t44 = t37 ^ 2 * t38;
t43 = pkin(3) + qJ(3);
t31 = qJD(2) * pkin(1);
t34 = sin(qJ(2));
t41 = qJD(1) * t34;
t18 = qJD(2) * pkin(2) - t43 * t41 + t31;
t40 = qJD(1) * t37;
t24 = t43 * t40;
t33 = sin(qJ(4));
t36 = cos(qJ(4));
t42 = t33 * t18 + t36 * t24;
t39 = qJD(1) * qJD(2);
t21 = t33 * t41 - t36 * t40;
t23 = qJD(3) + (-pkin(1) - pkin(2)) * t40;
t35 = cos(qJ(5));
t32 = sin(qJ(5));
t29 = qJD(2) + qJD(4);
t26 = -pkin(1) * t40 + qJD(3);
t25 = -qJ(3) * t41 + t31;
t22 = (t33 * t37 + t34 * t36) * qJD(1);
t19 = qJD(5) + t21;
t16 = t35 * t22 + t32 * t29;
t15 = t32 * t22 - t35 * t29;
t14 = -t22 * pkin(4) + t23;
t13 = -t36 * t18 + t33 * t24;
t12 = t29 * pkin(4) + t42;
t1 = [t45, 0, 0, t34 ^ 2 * t45, t34 * t38 * t37, t34 * t39, t37 * t39, qJD(2) ^ 2 / 0.2e1, 0, 0, t25 * qJD(2) - t26 * t40, (-qJ(3) * qJD(2) * t37 + t26 * t34) * qJD(1), qJ(3) * t44 - t25 * t41, qJ(3) ^ 2 * t44 / 0.2e1 + t25 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t21, t22 * t29, -t21 * t29, t29 ^ 2 / 0.2e1, -t13 * t29 + t23 * t21, t23 * t22 - t42 * t29, t16 ^ 2 / 0.2e1, -t16 * t15, t16 * t19, -t15 * t19, t19 ^ 2 / 0.2e1, (-t32 * t12 + t35 * t14) * t19 + t13 * t15, -(t35 * t12 + t32 * t14) * t19 + t13 * t16;];
T_reg = t1;
