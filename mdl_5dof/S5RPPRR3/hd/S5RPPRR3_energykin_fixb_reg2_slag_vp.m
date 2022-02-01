% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:31
% EndTime: 2022-01-23 09:14:31
% DurationCPUTime: 0.16s
% Computational Cost: add. (262->44), mult. (671->109), div. (0->0), fcn. (431->8), ass. (0->35)
t37 = qJD(1) ^ 2;
t30 = t37 / 0.2e1;
t44 = pkin(1) * t37;
t43 = cos(qJ(4));
t42 = cos(qJ(5));
t32 = sin(pkin(8));
t25 = (pkin(1) * t32 + qJ(3)) * qJD(1);
t33 = cos(pkin(9));
t28 = t33 * qJD(2);
t31 = sin(pkin(9));
t14 = t28 + (-pkin(6) * qJD(1) - t25) * t31;
t17 = t31 * qJD(2) + t33 * t25;
t40 = qJD(1) * t33;
t15 = pkin(6) * t40 + t17;
t36 = sin(qJ(4));
t6 = t36 * t14 + t43 * t15;
t41 = qJD(1) * t31;
t34 = cos(pkin(8));
t39 = -pkin(1) * t34 - pkin(2);
t5 = t43 * t14 - t36 * t15;
t19 = qJD(3) + (-pkin(3) * t33 + t39) * qJD(1);
t35 = sin(qJ(5));
t29 = qJD(4) + qJD(5);
t24 = t39 * qJD(1) + qJD(3);
t22 = (t43 * t31 + t33 * t36) * qJD(1);
t20 = t36 * t41 - t43 * t40;
t16 = -t31 * t25 + t28;
t10 = t20 * pkin(4) + t19;
t9 = -t35 * t20 + t42 * t22;
t7 = t42 * t20 + t35 * t22;
t4 = -t20 * pkin(7) + t6;
t3 = qJD(4) * pkin(4) - t22 * pkin(7) + t5;
t2 = t35 * t3 + t42 * t4;
t1 = t42 * t3 - t35 * t4;
t8 = [0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t34 * t44, -t32 * t44, 0, qJD(2) ^ 2 / 0.2e1 + (t32 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t37, t31 ^ 2 * t30, t31 * t37 * t33, 0, t33 ^ 2 * t30, 0, 0, -t24 * t40, t24 * t41, (-t16 * t31 + t17 * t33) * qJD(1), t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * qJD(4), t20 ^ 2 / 0.2e1, -t20 * qJD(4), qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) + t19 * t20, -t6 * qJD(4) + t19 * t22, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t29, t7 ^ 2 / 0.2e1, -t7 * t29, t29 ^ 2 / 0.2e1, t1 * t29 + t10 * t7, t10 * t9 - t2 * t29, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
