% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:42
% EndTime: 2019-12-05 15:19:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (213->36), mult. (575->99), div. (0->0), fcn. (464->12), ass. (0->38)
t22 = cos(pkin(5)) * qJD(1) + qJD(2);
t28 = sin(pkin(6));
t31 = cos(pkin(6));
t30 = cos(pkin(11));
t29 = sin(pkin(5));
t47 = qJD(1) * t29;
t42 = t30 * t47;
t51 = t22 * t28 + t31 * t42;
t34 = sin(qJ(3));
t36 = cos(qJ(3));
t27 = sin(pkin(11));
t43 = t27 * t47;
t11 = -t34 * t43 + t51 * t36;
t37 = qJD(3) ^ 2;
t50 = t37 / 0.2e1;
t12 = t51 * t34 + t36 * t43;
t10 = qJD(3) * pkin(8) + t12;
t14 = t31 * t22 - t28 * t42;
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t6 = t35 * t10 + t33 * t14;
t49 = cos(qJ(5));
t46 = qJD(3) * t33;
t45 = t35 * qJD(3);
t44 = qJD(3) * qJD(4);
t5 = -t33 * t10 + t35 * t14;
t38 = qJD(1) ^ 2;
t32 = sin(qJ(5));
t23 = -qJD(5) + t45;
t19 = t32 * qJD(4) + t49 * t46;
t17 = -t49 * qJD(4) + t32 * t46;
t9 = -qJD(3) * pkin(3) - t11;
t7 = (-pkin(4) * t35 - pkin(9) * t33 - pkin(3)) * qJD(3) - t11;
t4 = qJD(4) * pkin(9) + t6;
t3 = -qJD(4) * pkin(4) - t5;
t2 = t32 * t7 + t49 * t4;
t1 = -t32 * t4 + t49 * t7;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t38 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 ^ 2 / 0.2e1 + (t27 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1) * t38 * t29 ^ 2, 0, 0, 0, 0, 0, t50, t11 * qJD(3), -t12 * qJD(3), 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t33 ^ 2 * t50, t33 * t37 * t35, t33 * t44, t35 ^ 2 * t50, t35 * t44, qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) - t9 * t45, -t6 * qJD(4) + t9 * t46, (-t33 * t5 + t35 * t6) * qJD(3), t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t17, -t19 * t23, t17 ^ 2 / 0.2e1, t17 * t23, t23 ^ 2 / 0.2e1, -t1 * t23 + t3 * t17, t3 * t19 + t2 * t23, -t1 * t19 - t2 * t17, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
