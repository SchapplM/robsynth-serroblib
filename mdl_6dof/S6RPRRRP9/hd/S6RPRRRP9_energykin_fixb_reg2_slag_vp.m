% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:30
% EndTime: 2019-03-09 06:28:30
% DurationCPUTime: 0.14s
% Computational Cost: add. (483->55), mult. (961->115), div. (0->0), fcn. (579->6), ass. (0->43)
t48 = qJD(1) ^ 2;
t40 = t48 / 0.2e1;
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t26 = (pkin(3) * t45 - pkin(8) * t47 + qJ(2)) * qJD(1);
t35 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t27 = qJD(3) * pkin(8) + t45 * t35;
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t15 = t44 * t26 + t46 * t27;
t51 = qJD(1) * t47;
t29 = -t46 * qJD(3) + t44 * t51;
t11 = -t29 * pkin(9) + t15;
t43 = sin(qJ(5));
t53 = cos(qJ(5));
t14 = t46 * t26 - t44 * t27;
t31 = t44 * qJD(3) + t46 * t51;
t36 = t45 * qJD(1) + qJD(4);
t9 = t36 * pkin(4) - t31 * pkin(9) + t14;
t4 = t53 * t11 + t43 * t9;
t52 = t48 * qJ(2);
t50 = qJD(3) * t35;
t49 = qJD(1) * qJD(3);
t3 = -t43 * t11 + t53 * t9;
t28 = -qJD(3) * pkin(3) - t47 * t35;
t21 = t29 * pkin(4) + t28;
t42 = t47 ^ 2;
t41 = t45 ^ 2;
t38 = qJ(2) ^ 2 * t40;
t37 = -pkin(1) * qJD(1) + qJD(2);
t34 = qJD(5) + t36;
t32 = t34 ^ 2 / 0.2e1;
t20 = -t43 * t29 + t53 * t31;
t18 = t53 * t29 + t43 * t31;
t17 = t20 ^ 2 / 0.2e1;
t16 = t18 ^ 2 / 0.2e1;
t13 = t20 * t34;
t12 = t18 * t34;
t6 = t18 * pkin(5) + qJD(6) + t21;
t5 = t20 * t18;
t2 = -t18 * qJ(6) + t4;
t1 = t34 * pkin(5) - t20 * qJ(6) + t3;
t7 = [0, 0, 0, 0, 0, t40, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, t37 * qJD(1), t52, t38 + t37 ^ 2 / 0.2e1, t42 * t40, -t47 * t48 * t45, t47 * t49, t41 * t40, -t45 * t49, qJD(3) ^ 2 / 0.2e1, t45 * t52 + t47 * t50, -t45 * t50 + t47 * t52 (-t41 - t42) * t35 * qJD(1), t38 + (t41 / 0.2e1 + t42 / 0.2e1) * t35 ^ 2, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * t36, t29 ^ 2 / 0.2e1, -t29 * t36, t36 ^ 2 / 0.2e1, t14 * t36 + t28 * t29, -t15 * t36 + t28 * t31, -t14 * t31 - t15 * t29, t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t17, -t5, t13, t16, -t12, t32, t21 * t18 + t3 * t34, t21 * t20 - t4 * t34, -t4 * t18 - t3 * t20, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t17, -t5, t13, t16, -t12, t32, t1 * t34 + t6 * t18, -t2 * t34 + t6 * t20, -t1 * t20 - t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1;];
T_reg  = t7;
