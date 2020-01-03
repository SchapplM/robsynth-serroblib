% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:26
% EndTime: 2019-12-31 19:58:26
% DurationCPUTime: 0.17s
% Computational Cost: add. (333->48), mult. (867->107), div. (0->0), fcn. (559->6), ass. (0->43)
t42 = qJD(1) ^ 2;
t53 = t42 / 0.2e1;
t40 = sin(qJ(2));
t48 = qJD(1) * t40;
t50 = pkin(6) + qJ(3);
t30 = qJD(2) * pkin(2) - t50 * t48;
t41 = cos(qJ(2));
t47 = qJD(1) * t41;
t31 = t50 * t47;
t38 = sin(pkin(8));
t49 = cos(pkin(8));
t16 = t38 * t30 + t49 * t31;
t14 = qJD(2) * pkin(7) + t16;
t39 = sin(qJ(4));
t52 = cos(qJ(4));
t25 = t38 * t48 - t49 * t47;
t27 = (t38 * t41 + t49 * t40) * qJD(1);
t32 = qJD(3) + (-pkin(2) * t41 - pkin(1)) * qJD(1);
t9 = t25 * pkin(3) - t27 * pkin(7) + t32;
t4 = t52 * t14 + t39 * t9;
t51 = t41 * t42;
t46 = qJD(1) * qJD(2);
t3 = -t39 * t14 + t52 * t9;
t45 = t40 * t46;
t44 = t41 * t46;
t15 = t49 * t30 - t38 * t31;
t13 = -qJD(2) * pkin(3) - t15;
t37 = t41 ^ 2;
t36 = t40 ^ 2;
t35 = qJD(2) ^ 2 / 0.2e1;
t24 = qJD(4) + t25;
t22 = t24 ^ 2 / 0.2e1;
t21 = t39 * qJD(2) + t52 * t27;
t19 = -t52 * qJD(2) + t39 * t27;
t18 = t21 ^ 2 / 0.2e1;
t17 = t19 ^ 2 / 0.2e1;
t11 = t21 * t24;
t10 = t19 * t24;
t6 = t21 * t19;
t5 = t19 * pkin(4) + qJD(5) + t13;
t2 = -t19 * qJ(5) + t4;
t1 = t24 * pkin(4) - t21 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, t53, 0, 0, 0, 0, t36 * t53, t40 * t51, t45, t37 * t53, t44, t35, pkin(1) * t51 - pkin(6) * t45, -t42 * pkin(1) * t40 - pkin(6) * t44, (t36 + t37) * t42 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t37 / 0.2e1 + t36 / 0.2e1) * pkin(6) ^ 2) * t42, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * qJD(2), t25 ^ 2 / 0.2e1, -t25 * qJD(2), t35, t15 * qJD(2) + t32 * t25, -t16 * qJD(2) + t32 * t27, -t15 * t27 - t16 * t25, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t18, -t6, t11, t17, -t10, t22, t13 * t19 + t3 * t24, t13 * t21 - t4 * t24, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t18, -t6, t11, t17, -t10, t22, t1 * t24 + t5 * t19, -t2 * t24 + t5 * t21, -t1 * t21 - t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
