% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:47:41
% EndTime: 2020-01-03 11:47:41
% DurationCPUTime: 0.16s
% Computational Cost: add. (198->41), mult. (500->97), div. (0->0), fcn. (283->6), ass. (0->36)
t36 = qJD(1) ^ 2;
t30 = t36 / 0.2e1;
t31 = sin(pkin(8));
t23 = (pkin(1) * t31 + pkin(6)) * qJD(1);
t35 = cos(qJ(3));
t28 = t35 * qJD(2);
t34 = sin(qJ(3));
t10 = qJD(3) * pkin(3) + t28 + (-pkin(7) * qJD(1) - t23) * t34;
t15 = t34 * qJD(2) + t35 * t23;
t40 = qJD(1) * t35;
t11 = pkin(7) * t40 + t15;
t33 = sin(qJ(4));
t42 = cos(qJ(4));
t4 = t33 * t10 + t42 * t11;
t43 = pkin(1) * t36;
t41 = qJD(1) * t34;
t39 = qJD(1) * qJD(3);
t32 = cos(pkin(8));
t38 = -pkin(1) * t32 - pkin(2);
t3 = t42 * t10 - t33 * t11;
t21 = (-pkin(3) * t35 + t38) * qJD(1);
t29 = qJD(3) + qJD(4);
t26 = t29 ^ 2 / 0.2e1;
t24 = t38 * qJD(1);
t20 = (t33 * t35 + t42 * t34) * qJD(1);
t18 = t33 * t41 - t42 * t40;
t17 = t20 ^ 2 / 0.2e1;
t16 = t18 ^ 2 / 0.2e1;
t14 = -t34 * t23 + t28;
t13 = t20 * t29;
t12 = t18 * t29;
t6 = t20 * t18;
t5 = t18 * pkin(4) + qJD(5) + t21;
t2 = -t18 * qJ(5) + t4;
t1 = t29 * pkin(4) - t20 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t32 * t43, -t31 * t43, 0, qJD(2) ^ 2 / 0.2e1 + (t31 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t36, t34 ^ 2 * t30, t34 * t36 * t35, t34 * t39, t35 ^ 2 * t30, t35 * t39, qJD(3) ^ 2 / 0.2e1, t14 * qJD(3) - t24 * t40, -t15 * qJD(3) + t24 * t41, (-t14 * t34 + t15 * t35) * qJD(1), t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t17, -t6, t13, t16, -t12, t26, t21 * t18 + t3 * t29, t21 * t20 - t4 * t29, -t4 * t18 - t3 * t20, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t17, -t6, t13, t16, -t12, t26, t1 * t29 + t5 * t18, -t2 * t29 + t5 * t20, -t1 * t20 - t2 * t18, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
