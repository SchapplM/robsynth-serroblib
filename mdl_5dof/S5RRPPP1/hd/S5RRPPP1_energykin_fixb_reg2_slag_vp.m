% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:29
% EndTime: 2019-12-31 19:24:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (395->54), mult. (1159->111), div. (0->0), fcn. (803->6), ass. (0->45)
t38 = qJD(1) ^ 2;
t56 = t38 / 0.2e1;
t34 = sin(pkin(8));
t51 = cos(pkin(8));
t52 = cos(pkin(5));
t42 = t52 * t51;
t35 = sin(pkin(5));
t44 = t35 * t51;
t37 = cos(qJ(2));
t49 = qJD(1) * t37;
t36 = sin(qJ(2));
t50 = qJD(1) * t36;
t14 = -qJD(2) * t44 + t34 * t50 - t42 * t49;
t43 = t37 * t52;
t48 = qJD(2) * t35;
t16 = t34 * t48 - (-t34 * t43 - t51 * t36) * qJD(1);
t55 = t14 * t16;
t27 = -t52 * qJD(2) + t35 * t49;
t10 = t16 * t27;
t11 = t27 * t14;
t54 = t37 * t38;
t53 = pkin(3) + qJ(5);
t12 = t14 ^ 2 / 0.2e1;
t13 = t16 ^ 2 / 0.2e1;
t47 = qJD(1) * qJD(2);
t23 = pkin(7) * t49 + (qJD(1) * t43 + t48) * qJ(3);
t24 = qJD(2) * pkin(2) + (-t52 * qJ(3) - pkin(7)) * t50;
t25 = (-qJ(3) * t35 * t36 - pkin(2) * t37 - pkin(1)) * qJD(1);
t8 = t51 * t23 + (t24 * t52 + t25 * t35) * t34;
t46 = t36 * t47;
t45 = t37 * t47;
t9 = -t35 * t24 + t52 * t25 + qJD(3);
t6 = t27 * qJ(4) - t8;
t41 = -t16 * qJ(4) + t9;
t7 = -t34 * t23 + t24 * t42 + t25 * t44;
t40 = qJD(4) - t7;
t33 = t37 ^ 2;
t32 = t36 ^ 2;
t26 = t27 ^ 2 / 0.2e1;
t5 = t27 * pkin(3) + t40;
t4 = t14 * pkin(3) + t41;
t3 = -t14 * pkin(4) + qJD(5) - t6;
t2 = t53 * t14 + t41;
t1 = t16 * pkin(4) + t53 * t27 + t40;
t15 = [0, 0, 0, 0, 0, t56, 0, 0, 0, 0, t32 * t56, t36 * t54, t46, t33 * t56, t45, qJD(2) ^ 2 / 0.2e1, pkin(1) * t54 - pkin(7) * t46, -t38 * pkin(1) * t36 - pkin(7) * t45, (t32 + t33) * t38 * pkin(7), (pkin(1) ^ 2 / 0.2e1 + (t33 / 0.2e1 + t32 / 0.2e1) * pkin(7) ^ 2) * t38, t13, -t55, -t10, t12, t11, t26, t9 * t14 - t7 * t27, t9 * t16 + t8 * t27, -t8 * t14 - t7 * t16, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t26, t10, -t11, t13, -t55, t12, t6 * t14 + t5 * t16, -t4 * t14 - t5 * t27, -t4 * t16 + t6 * t27, t4 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t26, -t11, -t10, t12, t55, t13, t1 * t16 - t3 * t14, -t2 * t16 - t3 * t27, t1 * t27 + t2 * t14, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t15;
