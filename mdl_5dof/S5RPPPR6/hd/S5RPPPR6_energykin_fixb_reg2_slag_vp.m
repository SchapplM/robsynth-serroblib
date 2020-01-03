% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:49
% EndTime: 2019-12-31 17:47:50
% DurationCPUTime: 0.16s
% Computational Cost: add. (199->44), mult. (521->104), div. (0->0), fcn. (279->6), ass. (0->40)
t36 = qJD(1) ^ 2;
t47 = t36 / 0.2e1;
t32 = sin(pkin(7));
t40 = qJ(2) * qJD(1);
t18 = t32 * t40 + qJD(3);
t42 = qJD(1) * t32;
t15 = pkin(3) * t42 + t18;
t31 = sin(pkin(8));
t33 = cos(pkin(8));
t34 = cos(pkin(7));
t38 = -qJ(3) * t32 - pkin(1);
t9 = qJD(2) + ((-pkin(2) - qJ(4)) * t34 + t38) * qJD(1);
t6 = t31 * t15 + t33 * t9;
t46 = sin(qJ(5));
t30 = t34 ^ 2;
t45 = t30 * t36;
t44 = t31 * t34;
t43 = qJ(2) * t36;
t41 = qJD(1) * t34;
t19 = t32 * t36 * t34;
t23 = t45 / 0.2e1;
t39 = qJ(2) ^ 2 * t47;
t16 = pkin(3) * t41 + t34 * t40 + qJD(4);
t5 = t33 * t15 - t31 * t9;
t35 = cos(qJ(5));
t29 = t32 ^ 2;
t28 = -qJD(1) * pkin(1) + qJD(2);
t25 = t30 * t43;
t22 = t29 * t47;
t20 = t30 * t39;
t17 = t33 * t41 + qJD(5);
t14 = qJD(2) + (-pkin(2) * t34 + t38) * qJD(1);
t11 = t35 * t31 * t41 - t46 * t42;
t10 = (t32 * t35 + t46 * t44) * qJD(1);
t7 = (pkin(4) * t33 + pkin(6) * t31) * t41 + t16;
t4 = pkin(6) * t42 + t6;
t3 = -pkin(4) * t42 - t5;
t2 = t35 * t4 + t46 * t7;
t1 = t35 * t7 - t46 * t4;
t8 = [0, 0, 0, 0, 0, t47, 0, 0, 0, 0, t22, t19, 0, t23, 0, 0, -t28 * t41, t28 * t42, t29 * t43 + t25, t20 + t29 * t39 + t28 ^ 2 / 0.2e1, 0, 0, 0, t22, t19, t23, t18 * t42 + t25, t14 * t41, -t14 * t42, t14 ^ 2 / 0.2e1 + t20 + t18 ^ 2 / 0.2e1, t31 ^ 2 * t23, t31 * t33 * t45, -t31 * t19, t33 ^ 2 * t23, -t33 * t19, t22, (t16 * t33 * t34 + t32 * t5) * qJD(1), (-t16 * t44 - t32 * t6) * qJD(1), (t31 * t5 - t33 * t6) * t41, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t10, -t11 * t17, t10 ^ 2 / 0.2e1, t10 * t17, t17 ^ 2 / 0.2e1, t1 * t17 - t3 * t10, -t3 * t11 - t2 * t17, t1 * t11 + t2 * t10, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
