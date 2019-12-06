% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:31
% EndTime: 2019-12-05 17:51:31
% DurationCPUTime: 0.13s
% Computational Cost: add. (194->32), mult. (375->82), div. (0->0), fcn. (190->8), ass. (0->36)
t20 = qJD(1) + qJD(3);
t18 = t20 ^ 2;
t43 = t18 / 0.2e1;
t31 = qJD(1) ^ 2;
t42 = pkin(1) * t31;
t23 = sin(pkin(9));
t41 = t23 ^ 2 * t18;
t40 = t20 * t23;
t25 = cos(pkin(9));
t39 = t25 * t20;
t26 = cos(pkin(8));
t12 = (pkin(1) * t26 + pkin(2)) * qJD(1);
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t24 = sin(pkin(8));
t35 = pkin(1) * qJD(1) * t24;
t10 = t28 * t12 + t30 * t35;
t8 = t20 * qJ(4) + t10;
t4 = -t25 * qJD(2) + t23 * t8;
t38 = t4 ^ 2 / 0.2e1;
t27 = sin(qJ(5));
t37 = t27 * t40;
t29 = cos(qJ(5));
t36 = t29 * t40;
t34 = t41 / 0.2e1;
t9 = t30 * t12 - t28 * t35;
t33 = qJD(4) - t9;
t22 = t31 / 0.2e1;
t21 = qJD(2) ^ 2 / 0.2e1;
t13 = -qJD(5) + t39;
t7 = -t20 * pkin(3) + t33;
t6 = t23 * qJD(2) + t25 * t8;
t3 = (-pkin(4) * t25 - pkin(7) * t23 - pkin(3)) * t20 + t33;
t2 = t27 * t3 + t29 * t6;
t1 = -t27 * t6 + t29 * t3;
t5 = [0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t26 * t42, -t24 * t42, 0, t21 + (t24 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t31, 0, 0, 0, 0, 0, t43, t9 * t20, -t10 * t20, 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t21, t34, t23 * t18 * t25, 0, t25 ^ 2 * t43, 0, 0, -t7 * t39, t7 * t40, (t23 * t4 + t25 * t6) * t20, t6 ^ 2 / 0.2e1 + t38 + t7 ^ 2 / 0.2e1, t29 ^ 2 * t34, -t29 * t27 * t41, -t13 * t36, t27 ^ 2 * t34, t13 * t37, t13 ^ 2 / 0.2e1, -t1 * t13 + t4 * t37, t2 * t13 + t4 * t36, (-t1 * t29 - t2 * t27) * t40, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t38;];
T_reg = t5;
