% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:20
% EndTime: 2019-12-31 22:12:20
% DurationCPUTime: 0.19s
% Computational Cost: add. (535->52), mult. (1331->117), div. (0->0), fcn. (978->8), ass. (0->45)
t48 = sin(qJ(2));
t49 = cos(qJ(2));
t44 = sin(pkin(5));
t57 = qJD(1) * t44;
t52 = t49 * t57;
t56 = cos(pkin(5)) * qJD(1);
t54 = pkin(1) * t56;
t34 = pkin(7) * t52 + t48 * t54;
t42 = qJD(2) + t56;
t26 = pkin(8) * t42 + t34;
t28 = (-pkin(2) * t49 - pkin(8) * t48 - pkin(1)) * t57;
t47 = sin(qJ(3));
t60 = cos(qJ(3));
t16 = t60 * t26 + t47 * t28;
t37 = -qJD(3) + t52;
t12 = -pkin(9) * t37 + t16;
t46 = sin(qJ(4));
t59 = cos(qJ(4));
t53 = t48 * t57;
t33 = -pkin(7) * t53 + t49 * t54;
t25 = -pkin(2) * t42 - t33;
t30 = -t60 * t42 + t47 * t53;
t32 = t47 * t42 + t53 * t60;
t9 = pkin(3) * t30 - pkin(9) * t32 + t25;
t4 = t59 * t12 + t46 * t9;
t50 = qJD(1) ^ 2;
t58 = t44 ^ 2 * t50;
t55 = t49 * t58;
t51 = t58 / 0.2e1;
t3 = -t12 * t46 + t59 * t9;
t15 = -t47 * t26 + t28 * t60;
t11 = t37 * pkin(3) - t15;
t29 = qJD(4) + t30;
t27 = t29 ^ 2 / 0.2e1;
t21 = t32 * t59 - t46 * t37;
t19 = t46 * t32 + t37 * t59;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t14 = t21 * t29;
t13 = t19 * t29;
t6 = t21 * t19;
t5 = t19 * pkin(4) + qJD(5) + t11;
t2 = -qJ(5) * t19 + t4;
t1 = pkin(4) * t29 - qJ(5) * t21 + t3;
t7 = [0, 0, 0, 0, 0, t50 / 0.2e1, 0, 0, 0, 0, t48 ^ 2 * t51, t48 * t55, t42 * t53, t49 ^ 2 * t51, t42 * t52, t42 ^ 2 / 0.2e1, pkin(1) * t55 + t33 * t42, -pkin(1) * t48 * t58 - t34 * t42, (-t33 * t48 + t34 * t49) * t57, t34 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t51, t32 ^ 2 / 0.2e1, -t32 * t30, -t32 * t37, t30 ^ 2 / 0.2e1, t30 * t37, t37 ^ 2 / 0.2e1, -t15 * t37 + t25 * t30, t16 * t37 + t25 * t32, -t15 * t32 - t16 * t30, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t18, -t6, t14, t17, -t13, t27, t11 * t19 + t29 * t3, t11 * t21 - t29 * t4, -t19 * t4 - t21 * t3, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t18, -t6, t14, t17, -t13, t27, t1 * t29 + t19 * t5, -t2 * t29 + t21 * t5, -t1 * t21 - t19 * t2, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
