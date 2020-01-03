% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR14_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:15
% EndTime: 2019-12-31 18:35:15
% DurationCPUTime: 0.17s
% Computational Cost: add. (249->42), mult. (537->98), div. (0->0), fcn. (293->6), ass. (0->34)
t31 = sin(pkin(8));
t32 = cos(pkin(8));
t34 = sin(qJ(3));
t35 = cos(qJ(3));
t16 = (t31 * t35 + t32 * t34) * qJD(1);
t36 = qJD(1) ^ 2;
t27 = t36 / 0.2e1;
t42 = cos(qJ(5));
t21 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t38 = -qJ(4) * qJD(1) + t21;
t13 = qJD(3) * pkin(3) + t35 * t38;
t14 = t38 * t34;
t6 = t31 * t13 + t32 * t14;
t41 = t36 * qJ(2);
t40 = qJD(3) * t21;
t39 = qJD(1) * qJD(3);
t19 = qJD(4) + (pkin(3) * t34 + qJ(2)) * qJD(1);
t5 = t32 * t13 - t31 * t14;
t33 = sin(qJ(5));
t30 = t35 ^ 2;
t29 = t34 ^ 2;
t26 = qJD(3) ^ 2 / 0.2e1;
t24 = qJ(2) ^ 2 * t27;
t23 = -pkin(1) * qJD(1) + qJD(2);
t18 = (-t31 * t34 + t32 * t35) * qJD(1);
t15 = qJD(5) + t16;
t10 = t33 * qJD(3) + t42 * t18;
t8 = -t42 * qJD(3) + t33 * t18;
t7 = t16 * pkin(4) - t18 * pkin(7) + t19;
t4 = qJD(3) * pkin(7) + t6;
t3 = -qJD(3) * pkin(4) - t5;
t2 = t33 * t7 + t42 * t4;
t1 = -t33 * t4 + t42 * t7;
t9 = [0, 0, 0, 0, 0, t27, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, t23 * qJD(1), t41, t24 + t23 ^ 2 / 0.2e1, t30 * t27, -t35 * t36 * t34, t35 * t39, t29 * t27, -t34 * t39, t26, t34 * t41 + t35 * t40, -t34 * t40 + t35 * t41, (-t29 - t30) * t21 * qJD(1), t24 + (t29 / 0.2e1 + t30 / 0.2e1) * t21 ^ 2, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * qJD(3), t16 ^ 2 / 0.2e1, -t16 * qJD(3), t26, t5 * qJD(3) + t19 * t16, -t6 * qJD(3) + t19 * t18, -t6 * t16 - t5 * t18, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t15, t8 ^ 2 / 0.2e1, -t8 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t3 * t8, t3 * t10 - t2 * t15, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t9;
