% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:58:29
% EndTime: 2019-03-09 06:58:29
% DurationCPUTime: 0.16s
% Computational Cost: add. (682->61), mult. (1441->150), div. (0->0), fcn. (963->10), ass. (0->47)
t53 = qJD(1) ^ 2;
t44 = t53 / 0.2e1;
t61 = pkin(1) * t53;
t45 = sin(pkin(11));
t37 = (pkin(1) * t45 + pkin(7)) * qJD(1);
t52 = cos(qJ(3));
t42 = t52 * qJD(2);
t50 = sin(qJ(3));
t23 = qJD(3) * pkin(3) + t42 + (-pkin(8) * qJD(1) - t37) * t50;
t29 = t50 * qJD(2) + t52 * t37;
t57 = qJD(1) * t52;
t27 = pkin(8) * t57 + t29;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t15 = t49 * t23 + t51 * t27;
t43 = qJD(3) + qJD(4);
t10 = t43 * pkin(9) + t15;
t58 = qJD(1) * t50;
t32 = t49 * t58 - t51 * t57;
t34 = (t49 * t52 + t50 * t51) * qJD(1);
t46 = cos(pkin(11));
t55 = -pkin(1) * t46 - pkin(2);
t35 = (-pkin(3) * t52 + t55) * qJD(1);
t18 = t32 * pkin(4) - t34 * pkin(9) + t35;
t48 = sin(qJ(5));
t60 = cos(qJ(5));
t6 = t60 * t10 + t48 * t18;
t59 = cos(qJ(6));
t56 = qJD(1) * qJD(3);
t5 = -t48 * t10 + t60 * t18;
t14 = t51 * t23 - t49 * t27;
t9 = -t43 * pkin(4) - t14;
t31 = qJD(5) + t32;
t47 = sin(qJ(6));
t38 = t55 * qJD(1);
t30 = qJD(6) + t31;
t28 = -t50 * t37 + t42;
t26 = t60 * t34 + t48 * t43;
t24 = t48 * t34 - t60 * t43;
t13 = -t47 * t24 + t59 * t26;
t11 = t59 * t24 + t47 * t26;
t7 = t24 * pkin(5) + t9;
t4 = -t24 * pkin(10) + t6;
t3 = t31 * pkin(5) - t26 * pkin(10) + t5;
t2 = t47 * t3 + t59 * t4;
t1 = t59 * t3 - t47 * t4;
t8 = [0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t46 * t61, -t45 * t61, 0, qJD(2) ^ 2 / 0.2e1 + (t45 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t53, t50 ^ 2 * t44, t50 * t53 * t52, t50 * t56, t52 ^ 2 * t44, t52 * t56, qJD(3) ^ 2 / 0.2e1, t28 * qJD(3) - t38 * t57, -t29 * qJD(3) + t38 * t58 (-t28 * t50 + t29 * t52) * qJD(1), t29 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t34 ^ 2 / 0.2e1, -t34 * t32, t43 * t34, t32 ^ 2 / 0.2e1, -t43 * t32, t43 ^ 2 / 0.2e1, t14 * t43 + t35 * t32, -t15 * t43 + t35 * t34, -t14 * t34 - t15 * t32, t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, t26 * t31, t24 ^ 2 / 0.2e1, -t24 * t31, t31 ^ 2 / 0.2e1, t9 * t24 + t5 * t31, t9 * t26 - t6 * t31, -t6 * t24 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t30, t11 ^ 2 / 0.2e1, -t11 * t30, t30 ^ 2 / 0.2e1, t1 * t30 + t7 * t11, t7 * t13 - t2 * t30, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
