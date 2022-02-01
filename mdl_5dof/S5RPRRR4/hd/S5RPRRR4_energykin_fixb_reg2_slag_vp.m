% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:41
% EndTime: 2022-01-23 09:34:41
% DurationCPUTime: 0.10s
% Computational Cost: add. (173->26), mult. (347->72), div. (0->0), fcn. (172->8), ass. (0->30)
t16 = qJD(1) + qJD(3);
t15 = qJD(4) + t16;
t14 = t15 ^ 2;
t33 = t14 / 0.2e1;
t20 = cos(pkin(9));
t13 = (pkin(1) * t20 + pkin(2)) * qJD(1);
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t19 = sin(pkin(9));
t29 = pkin(1) * qJD(1) * t19;
t11 = t23 * t13 + t26 * t29;
t22 = sin(qJ(4));
t25 = cos(qJ(4));
t10 = t26 * t13 - t23 * t29;
t8 = t16 * pkin(3) + t10;
t6 = t25 * t11 + t22 * t8;
t27 = qJD(1) ^ 2;
t32 = pkin(1) * t27;
t5 = -t22 * t11 + t25 * t8;
t3 = -t15 * pkin(4) - t5;
t31 = t15 * t3;
t30 = qJD(5) * t15;
t24 = cos(qJ(5));
t21 = sin(qJ(5));
t18 = t27 / 0.2e1;
t17 = qJD(2) ^ 2 / 0.2e1;
t4 = t15 * pkin(8) + t6;
t2 = t21 * qJD(2) + t24 * t4;
t1 = t24 * qJD(2) - t21 * t4;
t7 = [0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t20 * t32, -t19 * t32, 0, t17 + (t19 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t27, 0, 0, 0, 0, 0, t16 ^ 2 / 0.2e1, t10 * t16, -t11 * t16, 0, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t17, 0, 0, 0, 0, 0, t33, t5 * t15, -t6 * t15, 0, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t17, t21 ^ 2 * t33, t21 * t14 * t24, t21 * t30, t24 ^ 2 * t33, t24 * t30, qJD(5) ^ 2 / 0.2e1, t1 * qJD(5) - t24 * t31, -t2 * qJD(5) + t21 * t31, (-t1 * t21 + t2 * t24) * t15, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
