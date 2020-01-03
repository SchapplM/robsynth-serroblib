% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:36
% EndTime: 2019-12-31 18:24:36
% DurationCPUTime: 0.15s
% Computational Cost: add. (158->41), mult. (386->96), div. (0->0), fcn. (179->6), ass. (0->38)
t31 = qJD(1) ^ 2;
t24 = t31 / 0.2e1;
t42 = -pkin(3) - pkin(7);
t41 = pkin(1) * t31;
t25 = sin(pkin(8));
t16 = (pkin(1) * t25 + pkin(6)) * qJD(1);
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t10 = t28 * qJD(2) + t30 * t16;
t40 = qJD(1) * t30;
t39 = t28 * qJD(1);
t38 = qJD(1) * qJD(3);
t26 = cos(pkin(8));
t37 = -pkin(1) * t26 - pkin(2);
t36 = t28 * t38;
t35 = t30 * t38;
t9 = t30 * qJD(2) - t28 * t16;
t7 = -qJD(3) * qJ(4) - t10;
t34 = qJD(4) - t9;
t33 = -qJ(4) * t28 + t37;
t29 = cos(qJ(5));
t27 = sin(qJ(5));
t23 = qJD(3) ^ 2 / 0.2e1;
t21 = t30 ^ 2 * t24;
t20 = t28 ^ 2 * t24;
t19 = qJD(5) + t39;
t18 = t28 * t31 * t30;
t17 = t37 * qJD(1);
t15 = t29 * qJD(3) - t27 * t40;
t13 = t27 * qJD(3) + t29 * t40;
t8 = (-pkin(3) * t30 + t33) * qJD(1);
t6 = -qJD(3) * pkin(3) + t34;
t5 = (t42 * t30 + t33) * qJD(1);
t4 = pkin(4) * t40 - t7;
t3 = pkin(4) * t39 + t42 * qJD(3) + t34;
t2 = t27 * t3 + t29 * t5;
t1 = -t27 * t5 + t29 * t3;
t11 = [0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t26 * t41, -t25 * t41, 0, qJD(2) ^ 2 / 0.2e1 + (t25 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t31, t20, t18, t36, t21, t35, t23, t9 * qJD(3) - t17 * t40, -t10 * qJD(3) + t17 * t39, (t10 * t30 - t28 * t9) * qJD(1), t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t23, -t36, -t35, t20, t18, t21, (t28 * t6 - t30 * t7) * qJD(1), t6 * qJD(3) + t8 * t40, -t7 * qJD(3) - t8 * t39, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t19, t13 ^ 2 / 0.2e1, -t13 * t19, t19 ^ 2 / 0.2e1, t1 * t19 + t4 * t13, t4 * t15 - t2 * t19, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg = t11;
