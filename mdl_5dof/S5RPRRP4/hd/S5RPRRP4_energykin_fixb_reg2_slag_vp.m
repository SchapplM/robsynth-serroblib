% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:23
% EndTime: 2020-01-03 11:50:23
% DurationCPUTime: 0.17s
% Computational Cost: add. (285->47), mult. (761->106), div. (0->0), fcn. (474->6), ass. (0->46)
t44 = qJD(1) ^ 2;
t57 = t44 / 0.2e1;
t38 = sin(pkin(8));
t39 = cos(pkin(8));
t23 = qJD(2) + (-pkin(2) * t39 - pkin(6) * t38 - pkin(1)) * qJD(1);
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t52 = qJ(2) * qJD(1);
t46 = t39 * t52;
t15 = t41 * t23 + t43 * t46;
t54 = qJD(1) * t38;
t50 = t41 * t54;
t11 = -pkin(7) * t50 + t15;
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t22 = t43 * t23;
t53 = t39 * qJD(1);
t32 = -qJD(3) + t53;
t7 = -t32 * pkin(3) + t22 + (-pkin(7) * t38 * t43 - qJ(2) * t39 * t41) * qJD(1);
t4 = t42 * t11 + t40 * t7;
t36 = t38 ^ 2;
t56 = t36 * t44;
t24 = pkin(3) * t50 + t38 * t52;
t55 = qJ(2) * t44;
t51 = t36 * t55;
t49 = t43 * t54;
t48 = t56 / 0.2e1;
t47 = qJ(2) ^ 2 * t57;
t3 = -t40 * t11 + t42 * t7;
t37 = t39 ^ 2;
t35 = -qJD(1) * pkin(1) + qJD(2);
t31 = t36 * t47;
t28 = -qJD(4) + t32;
t25 = t28 ^ 2 / 0.2e1;
t20 = -t40 * t50 + t42 * t49;
t18 = (-t40 * t43 - t41 * t42) * t54;
t17 = t20 ^ 2 / 0.2e1;
t16 = t18 ^ 2 / 0.2e1;
t14 = -t41 * t46 + t22;
t13 = t20 * t28;
t12 = t18 * t28;
t10 = -t18 * pkin(4) + qJD(5) + t24;
t8 = t20 * t18;
t2 = t18 * qJ(5) + t4;
t1 = -t28 * pkin(4) - t20 * qJ(5) + t3;
t5 = [0, 0, 0, 0, 0, t57, 0, 0, 0, 0, t48, t38 * t44 * t39, 0, t37 * t57, 0, 0, -t35 * t53, t35 * t54, (t36 + t37) * t55, t37 * t47 + t31 + t35 ^ 2 / 0.2e1, t43 ^ 2 * t48, -t43 * t41 * t56, -t32 * t49, t41 ^ 2 * t48, t32 * t50, t32 ^ 2 / 0.2e1, -t14 * t32 + t41 * t51, t15 * t32 + t43 * t51, (-t14 * t43 - t15 * t41) * t54, t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t31, t17, t8, -t13, t16, -t12, t25, -t24 * t18 - t3 * t28, t24 * t20 + t4 * t28, t4 * t18 - t3 * t20, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t17, t8, -t13, t16, -t12, t25, -t1 * t28 - t10 * t18, t10 * t20 + t2 * t28, -t1 * t20 + t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t5;
