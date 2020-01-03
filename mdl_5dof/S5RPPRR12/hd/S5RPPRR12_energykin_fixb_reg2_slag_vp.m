% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:20
% EndTime: 2019-12-31 18:07:20
% DurationCPUTime: 0.14s
% Computational Cost: add. (230->38), mult. (507->97), div. (0->0), fcn. (289->6), ass. (0->31)
t30 = sin(pkin(8));
t31 = cos(pkin(8));
t33 = sin(qJ(4));
t34 = cos(qJ(4));
t16 = (t30 * t34 + t31 * t33) * qJD(1);
t35 = qJD(1) ^ 2;
t28 = t35 / 0.2e1;
t39 = cos(qJ(5));
t21 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t37 = -pkin(6) * qJD(1) + t21;
t13 = t37 * t30;
t14 = t37 * t31;
t7 = t34 * t13 + t33 * t14;
t38 = qJD(1) * t30;
t23 = qJD(1) * qJ(2) + qJD(3);
t19 = pkin(3) * t38 + t23;
t6 = -t33 * t13 + t34 * t14;
t32 = sin(qJ(5));
t27 = t31 ^ 2;
t26 = t30 ^ 2;
t24 = -pkin(1) * qJD(1) + qJD(2);
t18 = (-t30 * t33 + t31 * t34) * qJD(1);
t15 = qJD(5) + t16;
t10 = t32 * qJD(4) + t39 * t18;
t8 = -t39 * qJD(4) + t32 * t18;
t5 = t16 * pkin(4) - t18 * pkin(7) + t19;
t4 = qJD(4) * pkin(7) + t7;
t3 = -qJD(4) * pkin(4) - t6;
t2 = t32 * t5 + t39 * t4;
t1 = -t32 * t4 + t39 * t5;
t9 = [0, 0, 0, 0, 0, t28, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, t24 * qJD(1), t35 * qJ(2), qJ(2) ^ 2 * t28 + t24 ^ 2 / 0.2e1, t27 * t28, -t31 * t35 * t30, 0, t26 * t28, 0, 0, t23 * t38, t23 * t31 * qJD(1), (-t26 - t27) * t21 * qJD(1), t23 ^ 2 / 0.2e1 + (t26 / 0.2e1 + t27 / 0.2e1) * t21 ^ 2, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * qJD(4), t16 ^ 2 / 0.2e1, -t16 * qJD(4), qJD(4) ^ 2 / 0.2e1, t6 * qJD(4) + t19 * t16, -t7 * qJD(4) + t19 * t18, -t7 * t16 - t6 * t18, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t15, t8 ^ 2 / 0.2e1, -t8 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t3 * t8, t3 * t10 - t2 * t15, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t9;
