% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:45:59
% EndTime: 2019-12-05 18:45:59
% DurationCPUTime: 0.18s
% Computational Cost: add. (388->48), mult. (1003->110), div. (0->0), fcn. (677->6), ass. (0->43)
t43 = qJD(1) ^ 2;
t53 = t43 / 0.2e1;
t52 = -pkin(7) - pkin(6);
t40 = sin(qJ(2));
t49 = qJD(1) * t40;
t29 = qJD(2) * pkin(2) + t52 * t49;
t42 = cos(qJ(2));
t48 = qJD(1) * t42;
t30 = t52 * t48;
t39 = sin(qJ(3));
t41 = cos(qJ(3));
t20 = t39 * t29 - t41 * t30;
t25 = t39 * t49 - t41 * t48;
t11 = -t25 * pkin(8) + t20;
t38 = sin(qJ(4));
t51 = cos(qJ(4));
t19 = t41 * t29 + t39 * t30;
t27 = (t39 * t42 + t40 * t41) * qJD(1);
t35 = qJD(2) + qJD(3);
t9 = t35 * pkin(3) - t27 * pkin(8) + t19;
t4 = t51 * t11 + t38 * t9;
t50 = t42 * t43;
t47 = qJD(1) * qJD(2);
t3 = -t38 * t11 + t51 * t9;
t46 = t40 * t47;
t45 = t42 * t47;
t31 = (-pkin(2) * t42 - pkin(1)) * qJD(1);
t21 = t25 * pkin(3) + t31;
t37 = t42 ^ 2;
t36 = t40 ^ 2;
t34 = qJD(4) + t35;
t33 = t34 ^ 2 / 0.2e1;
t18 = -t38 * t25 + t51 * t27;
t16 = t51 * t25 + t38 * t27;
t15 = t18 ^ 2 / 0.2e1;
t14 = t16 ^ 2 / 0.2e1;
t13 = t18 * t34;
t12 = t16 * t34;
t6 = t16 * pkin(4) + qJD(5) + t21;
t5 = t18 * t16;
t2 = -t16 * qJ(5) + t4;
t1 = t34 * pkin(4) - t18 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, t53, 0, 0, 0, 0, t36 * t53, t40 * t50, t46, t37 * t53, t45, qJD(2) ^ 2 / 0.2e1, pkin(1) * t50 - pkin(6) * t46, -t43 * pkin(1) * t40 - pkin(6) * t45, (t36 + t37) * t43 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t37 / 0.2e1 + t36 / 0.2e1) * pkin(6) ^ 2) * t43, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t35, t25 ^ 2 / 0.2e1, -t25 * t35, t35 ^ 2 / 0.2e1, t19 * t35 + t31 * t25, -t20 * t35 + t31 * t27, -t19 * t27 - t20 * t25, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t15, -t5, t13, t14, -t12, t33, t21 * t16 + t3 * t34, t21 * t18 - t4 * t34, -t4 * t16 - t3 * t18, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t15, -t5, t13, t14, -t12, t33, t1 * t34 + t6 * t16, t6 * t18 - t2 * t34, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg = t7;
