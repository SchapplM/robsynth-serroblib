% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:03
% EndTime: 2019-12-31 20:14:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (202->44), mult. (479->95), div. (0->0), fcn. (226->4), ass. (0->42)
t36 = qJD(1) ^ 2;
t50 = t36 / 0.2e1;
t49 = -pkin(2) - pkin(7);
t33 = sin(qJ(2));
t44 = t33 * qJD(1);
t42 = pkin(6) * t44 + qJD(3);
t10 = pkin(3) * t44 + t49 * qJD(2) + t42;
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t35 = cos(qJ(2));
t38 = -qJ(3) * t33 - pkin(1);
t8 = (t49 * t35 + t38) * qJD(1);
t4 = t32 * t10 + t34 * t8;
t45 = qJD(1) * t35;
t14 = t32 * qJD(2) + t34 * t45;
t16 = t34 * qJD(2) - t32 * t45;
t48 = t16 * t14;
t22 = qJD(4) + t44;
t47 = t22 * t14;
t46 = t35 * t36;
t18 = -pkin(6) * t45 - qJD(2) * qJ(3);
t43 = t14 ^ 2 / 0.2e1;
t41 = qJD(1) * qJD(2);
t11 = pkin(3) * t45 - t18;
t40 = t33 * t41;
t39 = t35 * t41;
t3 = t34 * t10 - t32 * t8;
t31 = t35 ^ 2;
t30 = t33 ^ 2;
t28 = qJD(2) ^ 2 / 0.2e1;
t24 = t31 * t50;
t23 = t30 * t50;
t21 = t33 * t46;
t19 = t22 ^ 2 / 0.2e1;
t17 = -qJD(2) * pkin(2) + t42;
t13 = t16 ^ 2 / 0.2e1;
t12 = (-pkin(2) * t35 + t38) * qJD(1);
t9 = t16 * t22;
t5 = t14 * pkin(4) - t16 * qJ(5) + t11;
t2 = t22 * qJ(5) + t4;
t1 = -t22 * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t50, 0, 0, 0, 0, t23, t21, t40, t24, t39, t28, pkin(1) * t46 - pkin(6) * t40, -t36 * pkin(1) * t33 - pkin(6) * t39, (t30 + t31) * t36 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t31 / 0.2e1 + t30 / 0.2e1) * pkin(6) ^ 2) * t36, t28, -t40, -t39, t23, t21, t24, (t17 * t33 - t18 * t35) * qJD(1), t17 * qJD(2) + t12 * t45, -t18 * qJD(2) - t12 * t44, t12 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t13, -t48, t9, t43, -t47, t19, t11 * t14 + t3 * t22, t11 * t16 - t4 * t22, -t4 * t14 - t3 * t16, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t13, t9, t48, t19, t47, t43, -t1 * t22 + t5 * t14, t1 * t16 - t2 * t14, -t5 * t16 + t2 * t22, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
