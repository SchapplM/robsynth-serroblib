% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR13_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:31:09
% EndTime: 2024-09-27 17:31:09
% DurationCPUTime: 0.05s
% Computational Cost: add. (483->36), mult. (678->105), div. (0->0), fcn. (401->10), ass. (0->44)
t28 = qJD(1) + qJD(2);
t37 = cos(qJ(2));
t46 = pkin(1) * qJD(1);
t41 = t37 * t46;
t22 = t28 * pkin(2) + t41;
t33 = sin(qJ(3));
t36 = cos(qJ(3));
t34 = sin(qJ(2));
t42 = t34 * t46;
t17 = t33 * t22 + t36 * t42;
t26 = qJD(3) + t28;
t29 = sin(pkin(5));
t48 = t26 * t29;
t11 = pkin(9) * t48 + t17;
t32 = sin(qJ(4));
t35 = cos(qJ(4));
t16 = t36 * t22 - t33 * t42;
t12 = t26 * pkin(3) + t16;
t30 = cos(pkin(5));
t50 = t12 * t30;
t6 = t35 * t11 + t32 * t50;
t51 = cos(qJ(5));
t25 = t26 ^ 2;
t27 = t29 ^ 2;
t49 = t25 * t27;
t47 = t26 * t35;
t45 = t12 * t26 * t27;
t44 = t32 * t48;
t43 = t29 * t47;
t40 = t49 / 0.2e1;
t23 = t30 * t26 + qJD(4);
t38 = qJD(1) ^ 2;
t31 = sin(qJ(5));
t21 = qJD(5) + t23;
t15 = (t31 * t35 + t51 * t32) * t48;
t13 = t31 * t44 - t51 * t43;
t9 = t35 * t50;
t7 = (-pkin(4) * t47 - t12) * t29;
t5 = -t32 * t11 + t9;
t4 = pkin(10) * t43 + t6;
t3 = t23 * pkin(4) + t9 + (-pkin(10) * t48 - t11) * t32;
t2 = t31 * t3 + t51 * t4;
t1 = t51 * t3 - t31 * t4;
t8 = [0, 0, 0, 0, 0, t38 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 ^ 2 / 0.2e1, t28 * t41, -t28 * t42, 0, (t34 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t38, 0, 0, 0, 0, 0, t25 / 0.2e1, t16 * t26, -t17 * t26, 0, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t32 ^ 2 * t40, t32 * t35 * t49, t23 * t44, t35 ^ 2 * t40, t23 * t43, t23 ^ 2 / 0.2e1, t5 * t23 + t35 * t45, -t6 * t23 - t32 * t45, (-t32 * t5 + t35 * t6) * t48, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t27 * t12 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t21, t13 ^ 2 / 0.2e1, -t13 * t21, t21 ^ 2 / 0.2e1, t1 * t21 + t7 * t13, t7 * t15 - t2 * t21, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
