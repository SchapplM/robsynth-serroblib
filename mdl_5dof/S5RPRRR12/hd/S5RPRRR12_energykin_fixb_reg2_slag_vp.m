% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:07
% EndTime: 2019-12-31 19:13:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (264->42), mult. (537->101), div. (0->0), fcn. (293->6), ass. (0->34)
t32 = sin(qJ(4));
t33 = sin(qJ(3));
t34 = cos(qJ(4));
t35 = cos(qJ(3));
t16 = (t32 * t35 + t33 * t34) * qJD(1);
t36 = qJD(1) ^ 2;
t27 = t36 / 0.2e1;
t42 = cos(qJ(5));
t21 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t38 = -pkin(7) * qJD(1) + t21;
t13 = qJD(3) * pkin(3) + t38 * t35;
t14 = t38 * t33;
t6 = t32 * t13 + t34 * t14;
t19 = (pkin(3) * t33 + qJ(2)) * qJD(1);
t41 = t36 * qJ(2);
t40 = qJD(3) * t21;
t39 = qJD(1) * qJD(3);
t5 = t34 * t13 - t32 * t14;
t31 = sin(qJ(5));
t30 = t35 ^ 2;
t29 = t33 ^ 2;
t26 = qJD(3) + qJD(4);
t25 = qJ(2) ^ 2 * t27;
t24 = -pkin(1) * qJD(1) + qJD(2);
t18 = (-t32 * t33 + t34 * t35) * qJD(1);
t15 = qJD(5) + t16;
t10 = t42 * t18 + t31 * t26;
t8 = t31 * t18 - t42 * t26;
t7 = t16 * pkin(4) - t18 * pkin(8) + t19;
t4 = t26 * pkin(8) + t6;
t3 = -t26 * pkin(4) - t5;
t2 = t31 * t7 + t42 * t4;
t1 = -t31 * t4 + t42 * t7;
t9 = [0, 0, 0, 0, 0, t27, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, t24 * qJD(1), t41, t25 + t24 ^ 2 / 0.2e1, t30 * t27, -t35 * t36 * t33, t35 * t39, t29 * t27, -t33 * t39, qJD(3) ^ 2 / 0.2e1, t33 * t41 + t35 * t40, -t33 * t40 + t35 * t41, (-t29 - t30) * t21 * qJD(1), t25 + (t29 / 0.2e1 + t30 / 0.2e1) * t21 ^ 2, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t26, t16 ^ 2 / 0.2e1, -t16 * t26, t26 ^ 2 / 0.2e1, t19 * t16 + t5 * t26, t19 * t18 - t6 * t26, -t6 * t16 - t5 * t18, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t15, t8 ^ 2 / 0.2e1, -t8 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t3 * t8, t3 * t10 - t2 * t15, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t9;
