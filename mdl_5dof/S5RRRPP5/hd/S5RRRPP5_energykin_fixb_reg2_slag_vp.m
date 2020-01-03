% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPP5
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:27
% EndTime: 2019-12-31 20:58:27
% DurationCPUTime: 0.18s
% Computational Cost: add. (227->43), mult. (583->93), div. (0->0), fcn. (335->4), ass. (0->39)
t29 = sin(qJ(3));
t30 = sin(qJ(2));
t31 = cos(qJ(2));
t44 = cos(qJ(3));
t16 = (t29 * t31 + t44 * t30) * qJD(1);
t21 = (pkin(2) * t31 + pkin(1)) * qJD(1);
t48 = -t16 * qJ(4) - t21;
t32 = qJD(1) ^ 2;
t47 = t32 / 0.2e1;
t46 = -pkin(3) - pkin(4);
t45 = -pkin(7) - pkin(6);
t39 = qJD(1) * t31;
t40 = qJD(1) * t30;
t14 = t29 * t40 - t44 * t39;
t9 = t16 * t14;
t25 = qJD(2) + qJD(3);
t10 = t16 * t25;
t43 = t25 * t14;
t42 = t31 * t32;
t19 = qJD(2) * pkin(2) + t45 * t40;
t20 = t45 * t39;
t8 = t29 * t19 - t44 * t20;
t12 = t14 ^ 2 / 0.2e1;
t23 = t25 ^ 2 / 0.2e1;
t38 = qJD(1) * qJD(2);
t6 = t25 * qJ(4) + t8;
t36 = t30 * t38;
t35 = t31 * t38;
t7 = t44 * t19 + t29 * t20;
t34 = qJD(4) - t7;
t28 = t31 ^ 2;
t27 = t30 ^ 2;
t13 = t16 ^ 2 / 0.2e1;
t5 = -t25 * pkin(3) + t34;
t4 = t14 * pkin(3) + t48;
t3 = t14 * qJ(5) + t6;
t2 = t46 * t14 + qJD(5) - t48;
t1 = -t16 * qJ(5) + t46 * t25 + t34;
t11 = [0, 0, 0, 0, 0, t47, 0, 0, 0, 0, t27 * t47, t30 * t42, t36, t28 * t47, t35, qJD(2) ^ 2 / 0.2e1, pkin(1) * t42 - pkin(6) * t36, -t32 * pkin(1) * t30 - pkin(6) * t35, (t27 + t28) * t32 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t28 / 0.2e1 + t27 / 0.2e1) * pkin(6) ^ 2) * t32, t13, -t9, t10, t12, -t43, t23, -t21 * t14 + t7 * t25, -t21 * t16 - t8 * t25, -t8 * t14 - t7 * t16, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t13, t10, t9, t23, t43, t12, t4 * t14 - t5 * t25, -t6 * t14 + t5 * t16, -t4 * t16 + t6 * t25, t6 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t13, t9, -t10, t12, -t43, t23, -t1 * t25 - t2 * t14, t2 * t16 + t3 * t25, -t1 * t16 + t3 * t14, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t11;
