% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:42
% EndTime: 2019-12-05 18:53:42
% DurationCPUTime: 0.13s
% Computational Cost: add. (188->34), mult. (432->100), div. (0->0), fcn. (267->8), ass. (0->42)
t26 = sin(qJ(3));
t21 = t26 ^ 2;
t48 = t21 / 0.2e1;
t28 = cos(qJ(3));
t23 = t28 ^ 2;
t47 = t23 / 0.2e1;
t46 = cos(qJ(4));
t45 = cos(qJ(5));
t20 = qJD(1) + qJD(2);
t44 = t20 * t28;
t29 = cos(qJ(2));
t43 = t20 * t29;
t30 = qJD(1) ^ 2;
t42 = t30 * pkin(1) ^ 2;
t41 = pkin(1) * qJD(1);
t25 = sin(qJ(4));
t27 = sin(qJ(2));
t37 = t27 * t41;
t32 = qJD(3) * pkin(2) - t26 * t37;
t34 = t28 * t37;
t6 = t25 * t34 - t46 * t32;
t40 = t6 ^ 2 / 0.2e1;
t39 = qJD(3) * t26;
t38 = qJD(3) * t28;
t36 = t29 * t41;
t35 = t42 / 0.2e1;
t33 = t20 * t37;
t10 = t25 * t26 * t20 - t46 * t44;
t24 = sin(qJ(5));
t22 = t27 ^ 2;
t19 = qJD(3) + qJD(4);
t18 = t20 ^ 2;
t16 = t29 ^ 2 * t35;
t14 = -pkin(2) * t44 - t36;
t12 = (t25 * t28 + t46 * t26) * t20;
t9 = qJD(5) + t10;
t8 = t25 * t32 + t46 * t34;
t5 = t45 * t12 + t24 * t19;
t3 = t24 * t12 - t45 * t19;
t2 = t24 * t14 + t45 * t8;
t1 = t45 * t14 - t24 * t8;
t4 = [0, 0, 0, 0, 0, t30 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 / 0.2e1, t20 * t36, -t33, 0, t22 * t35 + t16, t18 * t48, t26 * t18 * t28, t20 * t39, t18 * t47, t20 * t38, qJD(3) ^ 2 / 0.2e1, (-t27 * t39 + t28 * t43) * t41, (-t26 * t43 - t27 * t38) * t41, (t21 + t23) * t33, t16 + (t47 + t48) * t22 * t42, t12 ^ 2 / 0.2e1, -t12 * t10, t12 * t19, t10 ^ 2 / 0.2e1, -t10 * t19, t19 ^ 2 / 0.2e1, t14 * t10 - t6 * t19, t14 * t12 - t8 * t19, -t8 * t10 + t6 * t12, t8 ^ 2 / 0.2e1 + t40 + t14 ^ 2 / 0.2e1, t5 ^ 2 / 0.2e1, -t5 * t3, t5 * t9, t3 ^ 2 / 0.2e1, -t3 * t9, t9 ^ 2 / 0.2e1, t1 * t9 + t6 * t3, -t2 * t9 + t6 * t5, -t1 * t5 - t2 * t3, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t40;];
T_reg = t4;
