% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:51
% EndTime: 2019-12-05 16:06:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (146->34), mult. (393->80), div. (0->0), fcn. (228->4), ass. (0->31)
t27 = qJD(2) ^ 2;
t36 = t27 / 0.2e1;
t25 = sin(qJ(3));
t26 = cos(qJ(3));
t31 = qJD(2) * t26;
t16 = pkin(6) * t31 + t25 * qJD(1);
t11 = qJ(4) * t31 + t16;
t24 = sin(pkin(8));
t33 = cos(pkin(8));
t21 = t26 * qJD(1);
t32 = qJD(2) * t25;
t8 = qJD(3) * pkin(3) + t21 + (-pkin(6) - qJ(4)) * t32;
t5 = t33 * t11 + t24 * t8;
t12 = t24 * t32 - t33 * t31;
t14 = (t24 * t26 + t33 * t25) * qJD(2);
t35 = t14 * t12;
t34 = t26 * t27;
t30 = qJD(3) * t12;
t29 = t12 ^ 2 / 0.2e1;
t28 = qJD(2) * qJD(3);
t4 = -t24 * t11 + t33 * t8;
t17 = qJD(4) + (-pkin(3) * t26 - pkin(2)) * qJD(2);
t23 = qJD(1) ^ 2 / 0.2e1;
t22 = qJD(3) ^ 2 / 0.2e1;
t15 = -pkin(6) * t32 + t21;
t10 = t14 * qJD(3);
t9 = t14 ^ 2 / 0.2e1;
t3 = t12 * pkin(4) - t14 * qJ(5) + t17;
t2 = qJD(3) * qJ(5) + t5;
t1 = -qJD(3) * pkin(4) + qJD(5) - t4;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, t36, 0, 0, 0, t23, t25 ^ 2 * t36, t25 * t34, t25 * t28, t26 ^ 2 * t36, t26 * t28, t22, pkin(2) * t34 + t15 * qJD(3), -t27 * pkin(2) * t25 - t16 * qJD(3), (-t15 * t25 + t16 * t26) * qJD(2), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + pkin(2) ^ 2 * t36, t9, -t35, t10, t29, -t30, t22, t4 * qJD(3) + t17 * t12, -t5 * qJD(3) + t17 * t14, -t5 * t12 - t4 * t14, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t9, t10, t35, t22, t30, t29, -t1 * qJD(3) + t3 * t12, t1 * t14 - t2 * t12, t2 * qJD(3) - t3 * t14, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
