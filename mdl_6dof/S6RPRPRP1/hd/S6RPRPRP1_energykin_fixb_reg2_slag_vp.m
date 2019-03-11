% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:02
% EndTime: 2019-03-09 03:03:02
% DurationCPUTime: 0.16s
% Computational Cost: add. (462->57), mult. (1082->127), div. (0->0), fcn. (691->8), ass. (0->47)
t50 = qJD(1) ^ 2;
t43 = t50 / 0.2e1;
t59 = pkin(1) * t50;
t49 = cos(qJ(3));
t46 = cos(pkin(9));
t52 = -pkin(1) * t46 - pkin(2);
t31 = qJD(4) + (-pkin(3) * t49 + t52) * qJD(1);
t44 = sin(pkin(10));
t55 = qJD(1) * t49;
t48 = sin(qJ(3));
t56 = qJD(1) * t48;
t57 = cos(pkin(10));
t32 = t44 * t56 - t57 * t55;
t34 = (t44 * t49 + t57 * t48) * qJD(1);
t13 = t32 * pkin(4) - t34 * pkin(8) + t31;
t47 = sin(qJ(5));
t58 = cos(qJ(5));
t45 = sin(pkin(9));
t36 = (pkin(1) * t45 + pkin(7)) * qJD(1);
t41 = t49 * qJD(2);
t54 = qJ(4) * qJD(1);
t20 = qJD(3) * pkin(3) + t41 + (-t36 - t54) * t48;
t29 = t48 * qJD(2) + t49 * t36;
t26 = t49 * t54 + t29;
t10 = t44 * t20 + t57 * t26;
t8 = qJD(3) * pkin(8) + t10;
t4 = t47 * t13 + t58 * t8;
t53 = qJD(1) * qJD(3);
t3 = t58 * t13 - t47 * t8;
t9 = t57 * t20 - t44 * t26;
t7 = -qJD(3) * pkin(4) - t9;
t42 = qJD(3) ^ 2 / 0.2e1;
t37 = t52 * qJD(1);
t30 = qJD(5) + t32;
t28 = -t48 * t36 + t41;
t27 = t30 ^ 2 / 0.2e1;
t25 = t47 * qJD(3) + t58 * t34;
t23 = -t58 * qJD(3) + t47 * t34;
t22 = t25 ^ 2 / 0.2e1;
t21 = t23 ^ 2 / 0.2e1;
t16 = t25 * t30;
t15 = t23 * t30;
t14 = t25 * t23;
t5 = t23 * pkin(5) + qJD(6) + t7;
t2 = -t23 * qJ(6) + t4;
t1 = t30 * pkin(5) - t25 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t46 * t59, -t45 * t59, 0, qJD(2) ^ 2 / 0.2e1 + (t45 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t50, t48 ^ 2 * t43, t48 * t50 * t49, t48 * t53, t49 ^ 2 * t43, t49 * t53, t42, t28 * qJD(3) - t37 * t55, -t29 * qJD(3) + t37 * t56 (-t28 * t48 + t29 * t49) * qJD(1), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * qJD(3), t32 ^ 2 / 0.2e1, -t32 * qJD(3), t42, t9 * qJD(3) + t31 * t32, -t10 * qJD(3) + t31 * t34, -t10 * t32 - t9 * t34, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1, t22, -t14, t16, t21, -t15, t27, t7 * t23 + t3 * t30, t7 * t25 - t4 * t30, -t4 * t23 - t3 * t25, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t22, -t14, t16, t21, -t15, t27, t1 * t30 + t5 * t23, -t2 * t30 + t5 * t25, -t1 * t25 - t2 * t23, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
