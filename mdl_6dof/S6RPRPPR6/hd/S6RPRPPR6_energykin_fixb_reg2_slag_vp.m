% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:37
% EndTime: 2019-03-09 02:54:37
% DurationCPUTime: 0.18s
% Computational Cost: add. (597->58), mult. (1255->128), div. (0->0), fcn. (801->8), ass. (0->44)
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t48 = sin(qJ(3));
t49 = cos(qJ(3));
t29 = (t45 * t49 + t46 * t48) * qJD(1);
t50 = qJD(1) ^ 2;
t40 = t50 / 0.2e1;
t58 = cos(qJ(6));
t34 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t52 = -qJ(4) * qJD(1) + t34;
t26 = qJD(3) * pkin(3) + t52 * t49;
t27 = t52 * t48;
t17 = t45 * t26 + t46 * t27;
t13 = qJD(3) * qJ(5) + t17;
t31 = (-t45 * t48 + t46 * t49) * qJD(1);
t32 = qJD(4) + (pkin(3) * t48 + qJ(2)) * qJD(1);
t18 = t29 * pkin(4) - t31 * qJ(5) + t32;
t44 = sin(pkin(10));
t56 = cos(pkin(10));
t6 = t56 * t13 + t44 * t18;
t57 = t50 * qJ(2);
t55 = qJD(3) * t34;
t54 = t29 ^ 2 / 0.2e1;
t53 = qJD(1) * qJD(3);
t5 = -t44 * t13 + t56 * t18;
t16 = t46 * t26 - t45 * t27;
t12 = -qJD(3) * pkin(4) + qJD(5) - t16;
t47 = sin(qJ(6));
t43 = t49 ^ 2;
t42 = t48 ^ 2;
t39 = qJD(3) ^ 2 / 0.2e1;
t38 = qJ(2) ^ 2 * t40;
t36 = -pkin(1) * qJD(1) + qJD(2);
t28 = qJD(6) + t29;
t22 = t44 * qJD(3) + t56 * t31;
t20 = -t56 * qJD(3) + t44 * t31;
t10 = -t47 * t20 + t58 * t22;
t8 = t58 * t20 + t47 * t22;
t7 = t20 * pkin(5) + t12;
t4 = -t20 * pkin(8) + t6;
t3 = t29 * pkin(5) - t22 * pkin(8) + t5;
t2 = t47 * t3 + t58 * t4;
t1 = t58 * t3 - t47 * t4;
t9 = [0, 0, 0, 0, 0, t40, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, t36 * qJD(1), t57, t38 + t36 ^ 2 / 0.2e1, t43 * t40, -t49 * t50 * t48, t49 * t53, t42 * t40, -t48 * t53, t39, t48 * t57 + t49 * t55, -t48 * t55 + t49 * t57 (-t42 - t43) * t34 * qJD(1), t38 + (t42 / 0.2e1 + t43 / 0.2e1) * t34 ^ 2, t31 ^ 2 / 0.2e1, -t31 * t29, t31 * qJD(3), t54, -t29 * qJD(3), t39, t16 * qJD(3) + t32 * t29, -t17 * qJD(3) + t32 * t31, -t16 * t31 - t17 * t29, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t29, t20 ^ 2 / 0.2e1, -t20 * t29, t54, t12 * t20 + t5 * t29, t12 * t22 - t6 * t29, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t28, t8 ^ 2 / 0.2e1, -t8 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t8, t7 * t10 - t2 * t28, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t9;
