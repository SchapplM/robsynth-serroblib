% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR6
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:08:21
% EndTime: 2022-01-20 12:08:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (444->42), mult. (665->113), div. (0->0), fcn. (395->8), ass. (0->42)
t33 = sin(qJ(3));
t29 = t33 ^ 2;
t50 = t29 / 0.2e1;
t35 = cos(qJ(3));
t30 = t35 ^ 2;
t49 = t30 / 0.2e1;
t48 = cos(qJ(4));
t47 = cos(qJ(5));
t28 = qJD(1) + qJD(2);
t46 = t28 * t33;
t45 = t28 * t35;
t34 = sin(qJ(2));
t44 = pkin(1) * qJD(1);
t41 = t34 * t44;
t22 = t28 * pkin(7) + t41;
t39 = pkin(8) * t28 + t22;
t14 = qJD(3) * pkin(3) - t39 * t33;
t15 = t39 * t35;
t32 = sin(qJ(4));
t6 = t32 * t14 + t48 * t15;
t43 = qJD(3) * t33;
t42 = qJD(3) * t35;
t27 = qJD(3) + qJD(4);
t36 = cos(qJ(2));
t40 = t36 * t44;
t5 = t48 * t14 - t32 * t15;
t20 = -t40 + (-pkin(3) * t35 - pkin(2)) * t28;
t37 = qJD(1) ^ 2;
t31 = sin(qJ(5));
t26 = t28 ^ 2;
t25 = qJD(5) + t27;
t23 = -t28 * pkin(2) - t40;
t19 = (t32 * t35 + t48 * t33) * t28;
t17 = t32 * t46 - t48 * t45;
t10 = t17 * pkin(4) + t20;
t9 = -t31 * t17 + t47 * t19;
t7 = t47 * t17 + t31 * t19;
t4 = -t17 * pkin(9) + t6;
t3 = t27 * pkin(4) - t19 * pkin(9) + t5;
t2 = t31 * t3 + t47 * t4;
t1 = t47 * t3 - t31 * t4;
t8 = [0, 0, 0, 0, 0, t37 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 / 0.2e1, t28 * t40, -t28 * t41, 0, (t34 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t37, t26 * t50, t33 * t26 * t35, t28 * t43, t26 * t49, t28 * t42, qJD(3) ^ 2 / 0.2e1, -t22 * t43 - t23 * t45, -t22 * t42 + t23 * t46, (t29 + t30) * t28 * t22, t23 ^ 2 / 0.2e1 + (t49 + t50) * t22 ^ 2, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * t27, t17 ^ 2 / 0.2e1, -t17 * t27, t27 ^ 2 / 0.2e1, t20 * t17 + t5 * t27, t20 * t19 - t6 * t27, -t6 * t17 - t5 * t19, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t25, t7 ^ 2 / 0.2e1, -t7 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t10 * t7, t10 * t9 - t2 * t25, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
