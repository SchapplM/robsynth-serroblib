% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:40
% EndTime: 2019-12-31 20:53:40
% DurationCPUTime: 0.14s
% Computational Cost: add. (171->38), mult. (279->81), div. (0->0), fcn. (98->4), ass. (0->37)
t20 = sin(qJ(3));
t18 = t20 ^ 2;
t39 = t18 / 0.2e1;
t22 = cos(qJ(3));
t19 = t22 ^ 2;
t38 = t19 / 0.2e1;
t16 = qJD(1) + qJD(2);
t37 = t16 * t20;
t36 = t16 * t22;
t35 = -pkin(3) - qJ(5);
t34 = pkin(1) * qJD(1);
t21 = sin(qJ(2));
t29 = t21 * t34;
t9 = t16 * pkin(7) + t29;
t33 = t20 * t9 + qJD(4);
t32 = qJD(3) * t20;
t31 = qJD(3) * t22;
t30 = qJD(3) * qJ(4);
t23 = cos(qJ(2));
t28 = t23 * t34;
t27 = t16 * t31;
t26 = -qJ(4) * t20 - pkin(2);
t24 = qJD(1) ^ 2;
t17 = qJD(3) ^ 2 / 0.2e1;
t15 = t16 ^ 2;
t14 = t16 * t32;
t13 = t15 * t38;
t12 = t15 * t39;
t11 = t20 * t15 * t22;
t10 = -t16 * pkin(2) - t28;
t6 = -t22 * t9 - t30;
t5 = -qJD(3) * pkin(3) + t33;
t4 = -t28 + (-pkin(3) * t22 + t26) * t16;
t3 = t30 + qJD(5) + (pkin(4) * t16 + t9) * t22;
t2 = pkin(4) * t37 + t35 * qJD(3) + t33;
t1 = -t28 + (t35 * t22 + t26) * t16;
t7 = [0, 0, 0, 0, 0, t24 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 / 0.2e1, t16 * t28, -t16 * t29, 0, (t21 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t24, t12, t11, t14, t13, t27, t17, -t10 * t36 - t9 * t32, t10 * t37 - t9 * t31, (t18 + t19) * t9 * t16, t10 ^ 2 / 0.2e1 + (t38 + t39) * t9 ^ 2, t17, -t14, -t27, t12, t11, t13, (t20 * t5 - t22 * t6) * t16, t5 * qJD(3) + t4 * t36, -t6 * qJD(3) - t4 * t37, t4 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t17, -t27, t14, t13, -t11, t12, (t2 * t20 + t22 * t3) * t16, t3 * qJD(3) - t1 * t37, -t2 * qJD(3) - t1 * t36, t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t7;
