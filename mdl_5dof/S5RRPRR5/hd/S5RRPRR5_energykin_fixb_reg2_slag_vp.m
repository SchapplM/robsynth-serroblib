% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:48
% EndTime: 2022-01-20 11:02:48
% DurationCPUTime: 0.18s
% Computational Cost: add. (412->40), mult. (631->104), div. (0->0), fcn. (389->8), ass. (0->40)
t30 = sin(pkin(9));
t26 = t30 ^ 2;
t48 = t26 / 0.2e1;
t31 = cos(pkin(9));
t27 = t31 ^ 2;
t47 = t27 / 0.2e1;
t46 = cos(qJ(4));
t45 = cos(qJ(5));
t29 = qJD(1) + qJD(2);
t44 = t29 * t30;
t43 = t29 * t31;
t34 = sin(qJ(2));
t42 = pkin(1) * qJD(1);
t41 = t34 * t42;
t23 = t29 * qJ(3) + t41;
t39 = pkin(7) * t29 + t23;
t14 = t39 * t30;
t15 = t39 * t31;
t33 = sin(qJ(4));
t6 = -t33 * t14 + t46 * t15;
t35 = cos(qJ(2));
t40 = t35 * t42;
t5 = -t46 * t14 - t33 * t15;
t38 = qJD(3) - t40;
t17 = (-pkin(3) * t31 - pkin(2)) * t29 + t38;
t36 = qJD(1) ^ 2;
t32 = sin(qJ(5));
t28 = qJD(4) + qJD(5);
t25 = t29 ^ 2;
t21 = -t29 * pkin(2) + t38;
t20 = (t46 * t30 + t31 * t33) * t29;
t18 = t33 * t44 - t46 * t43;
t10 = t18 * pkin(4) + t17;
t9 = -t32 * t18 + t45 * t20;
t7 = t45 * t18 + t32 * t20;
t4 = -t18 * pkin(8) + t6;
t3 = qJD(4) * pkin(4) - t20 * pkin(8) + t5;
t2 = t32 * t3 + t45 * t4;
t1 = t45 * t3 - t32 * t4;
t8 = [0, 0, 0, 0, 0, t36 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 / 0.2e1, t29 * t40, -t29 * t41, 0, (t34 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t36, t25 * t48, t30 * t25 * t31, 0, t25 * t47, 0, 0, -t21 * t43, t21 * t44, (t26 + t27) * t29 * t23, t21 ^ 2 / 0.2e1 + (t47 + t48) * t23 ^ 2, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * qJD(4), t18 ^ 2 / 0.2e1, -t18 * qJD(4), qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) + t17 * t18, -t6 * qJD(4) + t17 * t20, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t28, t7 ^ 2 / 0.2e1, -t7 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t10 * t7, t10 * t9 - t2 * t28, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
