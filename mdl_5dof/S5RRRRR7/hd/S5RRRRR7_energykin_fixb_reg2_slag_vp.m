% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:25
% EndTime: 2019-12-31 22:22:26
% DurationCPUTime: 0.17s
% Computational Cost: add. (515->52), mult. (1253->130), div. (0->0), fcn. (870->8), ass. (0->44)
t43 = qJD(1) ^ 2;
t55 = t43 / 0.2e1;
t54 = -pkin(7) - pkin(6);
t41 = sin(qJ(2));
t49 = qJD(1) * t41;
t29 = qJD(2) * pkin(2) + t54 * t49;
t42 = cos(qJ(2));
t48 = qJD(1) * t42;
t30 = t54 * t48;
t40 = sin(qJ(3));
t53 = cos(qJ(3));
t20 = t40 * t29 - t53 * t30;
t25 = t40 * t49 - t53 * t48;
t11 = -t25 * pkin(8) + t20;
t39 = sin(qJ(4));
t52 = cos(qJ(4));
t19 = t53 * t29 + t40 * t30;
t27 = (t40 * t42 + t53 * t41) * qJD(1);
t35 = qJD(2) + qJD(3);
t9 = t35 * pkin(3) - t27 * pkin(8) + t19;
t6 = t52 * t11 + t39 * t9;
t51 = cos(qJ(5));
t50 = t42 * t43;
t47 = qJD(1) * qJD(2);
t46 = t41 * t47;
t45 = t42 * t47;
t16 = t52 * t25 + t39 * t27;
t31 = (-pkin(2) * t42 - pkin(1)) * qJD(1);
t5 = -t39 * t11 + t52 * t9;
t21 = t25 * pkin(3) + t31;
t38 = sin(qJ(5));
t37 = t42 ^ 2;
t36 = t41 ^ 2;
t34 = qJD(4) + t35;
t18 = -t39 * t25 + t52 * t27;
t15 = qJD(5) + t16;
t14 = t51 * t18 + t38 * t34;
t12 = t38 * t18 - t51 * t34;
t7 = t16 * pkin(4) - t18 * pkin(9) + t21;
t4 = t34 * pkin(9) + t6;
t3 = -t34 * pkin(4) - t5;
t2 = t38 * t7 + t51 * t4;
t1 = -t38 * t4 + t51 * t7;
t8 = [0, 0, 0, 0, 0, t55, 0, 0, 0, 0, t36 * t55, t41 * t50, t46, t37 * t55, t45, qJD(2) ^ 2 / 0.2e1, pkin(1) * t50 - pkin(6) * t46, -t43 * pkin(1) * t41 - pkin(6) * t45, (t36 + t37) * t43 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t37 / 0.2e1 + t36 / 0.2e1) * pkin(6) ^ 2) * t43, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t35, t25 ^ 2 / 0.2e1, -t25 * t35, t35 ^ 2 / 0.2e1, t19 * t35 + t31 * t25, -t20 * t35 + t31 * t27, -t19 * t27 - t20 * t25, t20 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t34, t16 ^ 2 / 0.2e1, -t16 * t34, t34 ^ 2 / 0.2e1, t21 * t16 + t5 * t34, t21 * t18 - t6 * t34, -t6 * t16 - t5 * t18, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t14 ^ 2 / 0.2e1, -t14 * t12, t14 * t15, t12 ^ 2 / 0.2e1, -t12 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t3 * t12, t3 * t14 - t2 * t15, -t1 * t14 - t2 * t12, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
