% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:56
% EndTime: 2019-12-05 18:42:57
% DurationCPUTime: 0.15s
% Computational Cost: add. (434->42), mult. (665->110), div. (0->0), fcn. (395->8), ass. (0->42)
t33 = sin(qJ(3));
t29 = t33 ^ 2;
t50 = t29 / 0.2e1;
t35 = cos(qJ(3));
t30 = t35 ^ 2;
t49 = t30 / 0.2e1;
t48 = cos(qJ(5));
t27 = qJD(1) + qJD(2);
t47 = t27 * t33;
t46 = t27 * t35;
t34 = sin(qJ(2));
t45 = pkin(1) * qJD(1);
t41 = t34 * t45;
t22 = t27 * pkin(7) + t41;
t39 = qJ(4) * t27 + t22;
t14 = qJD(3) * pkin(3) - t39 * t33;
t16 = t39 * t35;
t31 = sin(pkin(9));
t44 = cos(pkin(9));
t6 = t31 * t14 + t44 * t16;
t43 = qJD(3) * t33;
t42 = qJD(3) * t35;
t36 = cos(qJ(2));
t40 = t36 * t45;
t5 = t44 * t14 - t31 * t16;
t17 = -t40 + qJD(4) + (-pkin(3) * t35 - pkin(2)) * t27;
t37 = qJD(1) ^ 2;
t32 = sin(qJ(5));
t28 = qJD(3) ^ 2 / 0.2e1;
t26 = qJD(3) + qJD(5);
t25 = t27 ^ 2;
t23 = -t27 * pkin(2) - t40;
t20 = (t31 * t35 + t44 * t33) * t27;
t18 = t31 * t47 - t44 * t46;
t10 = t18 * pkin(4) + t17;
t9 = -t32 * t18 + t48 * t20;
t7 = t48 * t18 + t32 * t20;
t4 = -t18 * pkin(8) + t6;
t3 = qJD(3) * pkin(4) - t20 * pkin(8) + t5;
t2 = t32 * t3 + t48 * t4;
t1 = t48 * t3 - t32 * t4;
t8 = [0, 0, 0, 0, 0, t37 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 / 0.2e1, t27 * t40, -t27 * t41, 0, (t34 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t37, t25 * t50, t33 * t25 * t35, t27 * t43, t25 * t49, t27 * t42, t28, -t22 * t43 - t23 * t46, -t22 * t42 + t23 * t47, (t29 + t30) * t27 * t22, t23 ^ 2 / 0.2e1 + (t49 + t50) * t22 ^ 2, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * qJD(3), t18 ^ 2 / 0.2e1, -t18 * qJD(3), t28, t5 * qJD(3) + t17 * t18, -t6 * qJD(3) + t17 * t20, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t26, t7 ^ 2 / 0.2e1, -t7 * t26, t26 ^ 2 / 0.2e1, t1 * t26 + t10 * t7, t10 * t9 - t2 * t26, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
