% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR12_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:19:57
% EndTime: 2019-03-09 04:19:57
% DurationCPUTime: 0.14s
% Computational Cost: add. (373->58), mult. (714->120), div. (0->0), fcn. (357->6), ass. (0->45)
t45 = qJD(1) ^ 2;
t36 = t45 / 0.2e1;
t56 = cos(qJ(5));
t55 = cos(qJ(6));
t27 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t44 = cos(qJ(3));
t13 = qJD(4) + (pkin(4) * qJD(1) - t27) * t44 + (-pkin(3) - pkin(8)) * qJD(3);
t43 = sin(qJ(3));
t52 = qJD(1) * t43;
t54 = pkin(3) * t52 + qJD(1) * qJ(2);
t15 = (pkin(8) * t43 - qJ(4) * t44) * qJD(1) + t54;
t42 = sin(qJ(5));
t6 = t42 * t13 + t56 * t15;
t20 = -qJD(3) * qJ(4) - t43 * t27;
t53 = t45 * qJ(2);
t51 = qJD(3) * t27;
t50 = t44 * qJD(1);
t49 = qJD(1) * qJD(3);
t48 = t44 * t45 * t43;
t47 = t43 * t49;
t46 = t44 * t49;
t5 = t56 * t13 - t42 * t15;
t29 = qJD(5) + t50;
t16 = -pkin(4) * t52 - t20;
t41 = sin(qJ(6));
t40 = t44 ^ 2;
t39 = t43 ^ 2;
t35 = qJD(3) ^ 2 / 0.2e1;
t34 = qJ(2) ^ 2 * t36;
t33 = -pkin(1) * qJD(1) + qJD(2);
t31 = t40 * t36;
t30 = t39 * t36;
t26 = qJD(6) + t29;
t23 = t56 * qJD(3) + t42 * t52;
t21 = t42 * qJD(3) - t56 * t52;
t19 = -qJ(4) * t50 + t54;
t17 = -qJD(3) * pkin(3) - t44 * t27 + qJD(4);
t10 = -t41 * t21 + t55 * t23;
t8 = t55 * t21 + t41 * t23;
t7 = t21 * pkin(5) + t16;
t4 = -t21 * pkin(9) + t6;
t3 = t29 * pkin(5) - t23 * pkin(9) + t5;
t2 = t41 * t3 + t55 * t4;
t1 = t55 * t3 - t41 * t4;
t9 = [0, 0, 0, 0, 0, t36, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, t33 * qJD(1), t53, t34 + t33 ^ 2 / 0.2e1, t31, -t48, t46, t30, -t47, t35, t43 * t53 + t44 * t51, -t43 * t51 + t44 * t53 (-t39 - t40) * t27 * qJD(1), t34 + (t39 / 0.2e1 + t40 / 0.2e1) * t27 ^ 2, t35, -t46, t47, t31, -t48, t30 (t17 * t44 + t20 * t43) * qJD(1), t17 * qJD(3) - t19 * t52, -t20 * qJD(3) - t19 * t50, t19 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * t29, t21 ^ 2 / 0.2e1, -t21 * t29, t29 ^ 2 / 0.2e1, t16 * t21 + t5 * t29, t16 * t23 - t6 * t29, -t6 * t21 - t5 * t23, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t26, t8 ^ 2 / 0.2e1, -t8 * t26, t26 ^ 2 / 0.2e1, t1 * t26 + t7 * t8, t7 * t10 - t2 * t26, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t9;
