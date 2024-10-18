% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:11
% EndTime: 2024-09-28 18:08:11
% DurationCPUTime: 0.03s
% Computational Cost: add. (257->29), mult. (517->84), div. (0->0), fcn. (376->12), ass. (0->41)
t35 = qJD(1) ^ 2;
t46 = t35 / 0.2e1;
t34 = cos(qJ(2));
t24 = sin(pkin(5));
t43 = qJD(1) * t24;
t15 = qJD(2) * pkin(2) + t34 * t43;
t29 = sin(qJ(3));
t33 = cos(qJ(3));
t30 = sin(qJ(2));
t39 = t30 * t43;
t12 = t33 * t15 - t29 * t39;
t22 = qJD(2) + qJD(3);
t10 = t22 * pkin(3) + t12;
t13 = t29 * t15 + t33 * t39;
t28 = sin(qJ(4));
t32 = cos(qJ(4));
t7 = t28 * t10 + t32 * t13;
t19 = qJD(4) + t22;
t18 = t19 ^ 2;
t23 = sin(pkin(6));
t45 = t18 * t23 ^ 2;
t44 = t19 * t23;
t26 = cos(pkin(5));
t42 = qJD(1) * t26;
t27 = sin(qJ(5));
t41 = t27 * t44;
t31 = cos(qJ(5));
t40 = t31 * t44;
t38 = t45 / 0.2e1;
t6 = t32 * t10 - t28 * t13;
t37 = qJD(2) * t43;
t25 = cos(pkin(6));
t5 = t19 * pkin(4) + t6;
t36 = t23 * t42 + t25 * t5;
t17 = t26 ^ 2 * t46;
t16 = t25 * t19 + qJD(5);
t4 = pkin(10) * t44 + t7;
t3 = -t23 * t5 + t25 * t42;
t2 = t36 * t27 + t31 * t4;
t1 = -t27 * t4 + t36 * t31;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, qJD(2) ^ 2 / 0.2e1, t34 * t37, -t30 * t37, 0, t17 + (t30 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1) * t35 * t24 ^ 2, 0, 0, 0, 0, 0, t22 ^ 2 / 0.2e1, t12 * t22, -t13 * t22, 0, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t17, 0, 0, 0, 0, 0, t18 / 0.2e1, t6 * t19, -t7 * t19, 0, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t17, t27 ^ 2 * t38, t27 * t31 * t45, t16 * t41, t31 ^ 2 * t38, t16 * t40, t16 ^ 2 / 0.2e1, t1 * t16 - t3 * t40, -t2 * t16 + t3 * t41, (-t1 * t27 + t2 * t31) * t44, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
