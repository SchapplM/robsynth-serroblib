% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:56
% EndTime: 2019-12-31 20:49:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (290->36), mult. (458->91), div. (0->0), fcn. (239->6), ass. (0->40)
t26 = sin(qJ(3));
t23 = t26 ^ 2;
t45 = t23 / 0.2e1;
t28 = cos(qJ(3));
t24 = t28 ^ 2;
t44 = t24 / 0.2e1;
t25 = sin(pkin(8));
t39 = cos(pkin(8));
t21 = qJD(1) + qJD(2);
t27 = sin(qJ(2));
t40 = pkin(1) * qJD(1);
t34 = t27 * t40;
t17 = t21 * pkin(7) + t34;
t32 = qJ(4) * t21 + t17;
t8 = qJD(3) * pkin(3) - t32 * t26;
t9 = t32 * t28;
t5 = t25 * t8 + t39 * t9;
t41 = t21 * t28;
t42 = t21 * t26;
t13 = t25 * t42 - t39 * t41;
t15 = (t25 * t28 + t39 * t26) * t21;
t43 = t15 * t13;
t38 = qJD(3) * t13;
t37 = qJD(3) * t26;
t36 = qJD(3) * t28;
t35 = t13 ^ 2 / 0.2e1;
t29 = cos(qJ(2));
t33 = t29 * t40;
t4 = -t25 * t9 + t39 * t8;
t12 = -t33 + qJD(4) + (-pkin(3) * t28 - pkin(2)) * t21;
t30 = qJD(1) ^ 2;
t22 = qJD(3) ^ 2 / 0.2e1;
t20 = t21 ^ 2;
t18 = -t21 * pkin(2) - t33;
t11 = t15 * qJD(3);
t10 = t15 ^ 2 / 0.2e1;
t3 = qJD(3) * qJ(5) + t5;
t2 = -qJD(3) * pkin(4) + qJD(5) - t4;
t1 = t13 * pkin(4) - t15 * qJ(5) + t12;
t6 = [0, 0, 0, 0, 0, t30 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 / 0.2e1, t21 * t33, -t21 * t34, 0, (t27 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t30, t20 * t45, t26 * t20 * t28, t21 * t37, t20 * t44, t21 * t36, t22, -t17 * t37 - t18 * t41, -t17 * t36 + t18 * t42, (t23 + t24) * t21 * t17, t18 ^ 2 / 0.2e1 + (t44 + t45) * t17 ^ 2, t10, -t43, t11, t35, -t38, t22, t4 * qJD(3) + t12 * t13, -t5 * qJD(3) + t12 * t15, -t5 * t13 - t4 * t15, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10, t11, t43, t22, t38, t35, -t2 * qJD(3) + t1 * t13, -t3 * t13 + t2 * t15, t3 * qJD(3) - t1 * t15, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t6;
