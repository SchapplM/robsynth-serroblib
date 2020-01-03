% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR16_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:24
% EndTime: 2019-12-31 20:47:25
% DurationCPUTime: 0.17s
% Computational Cost: add. (406->55), mult. (1035->123), div. (0->0), fcn. (701->8), ass. (0->49)
t58 = -pkin(2) - pkin(8);
t44 = cos(qJ(2));
t42 = sin(qJ(2));
t48 = -qJ(3) * t42 - pkin(1);
t38 = sin(pkin(5));
t55 = qJD(1) * t38;
t15 = (t58 * t44 + t48) * t55;
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t50 = t42 * t55;
t32 = pkin(7) * t50;
t39 = cos(pkin(5));
t54 = t39 * qJD(1);
t36 = qJD(2) + t54;
t9 = qJD(3) + t32 + t58 * t36 + (-pkin(1) * t39 * t44 + pkin(3) * t38 * t42) * qJD(1);
t6 = t43 * t15 + t41 * t9;
t57 = cos(qJ(5));
t45 = qJD(1) ^ 2;
t56 = t38 ^ 2 * t45;
t51 = t44 * t55;
t52 = pkin(1) * t54;
t24 = pkin(7) * t51 + t42 * t52;
t53 = t44 * t56;
t17 = -t36 * qJ(3) - t24;
t49 = t56 / 0.2e1;
t14 = pkin(3) * t51 - t17;
t47 = t36 * t50;
t46 = t36 * t51;
t5 = -t41 * t15 + t43 * t9;
t23 = t44 * t52 - t32;
t20 = t41 * t36 + t43 * t51;
t40 = sin(qJ(5));
t30 = t44 ^ 2 * t49;
t29 = t42 ^ 2 * t49;
t28 = t36 ^ 2 / 0.2e1;
t27 = qJD(4) + t50;
t26 = t42 * t53;
t22 = t43 * t36 - t41 * t51;
t19 = qJD(5) + t20;
t18 = (-pkin(2) * t44 + t48) * t55;
t16 = -t36 * pkin(2) + qJD(3) - t23;
t12 = t57 * t22 + t40 * t27;
t10 = t40 * t22 - t57 * t27;
t7 = t20 * pkin(4) - t22 * pkin(9) + t14;
t4 = t27 * pkin(9) + t6;
t3 = -t27 * pkin(4) - t5;
t2 = t57 * t4 + t40 * t7;
t1 = -t40 * t4 + t57 * t7;
t8 = [0, 0, 0, 0, 0, t45 / 0.2e1, 0, 0, 0, 0, t29, t26, t47, t30, t46, t28, pkin(1) * t53 + t23 * t36, -pkin(1) * t42 * t56 - t24 * t36, (-t23 * t42 + t24 * t44) * t55, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t49, t28, -t47, -t46, t29, t26, t30, (t16 * t42 - t17 * t44) * t55, t16 * t36 + t18 * t51, -t17 * t36 - t18 * t50, t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t27, t20 ^ 2 / 0.2e1, -t20 * t27, t27 ^ 2 / 0.2e1, t14 * t20 + t5 * t27, t14 * t22 - t6 * t27, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t12 ^ 2 / 0.2e1, -t12 * t10, t12 * t19, t10 ^ 2 / 0.2e1, -t10 * t19, t19 ^ 2 / 0.2e1, t1 * t19 + t3 * t10, t3 * t12 - t2 * t19, -t1 * t12 - t2 * t10, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
