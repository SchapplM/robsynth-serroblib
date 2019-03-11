% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:38:56
% EndTime: 2019-03-09 03:38:56
% DurationCPUTime: 0.16s
% Computational Cost: add. (629->61), mult. (1441->147), div. (0->0), fcn. (963->10), ass. (0->48)
t53 = qJD(1) ^ 2;
t44 = t53 / 0.2e1;
t62 = pkin(1) * t53;
t46 = sin(pkin(10));
t37 = (pkin(1) * t46 + pkin(7)) * qJD(1);
t52 = cos(qJ(3));
t42 = t52 * qJD(2);
t51 = sin(qJ(3));
t57 = qJ(4) * qJD(1);
t23 = qJD(3) * pkin(3) + t42 + (-t37 - t57) * t51;
t30 = t51 * qJD(2) + t52 * t37;
t27 = t52 * t57 + t30;
t45 = sin(pkin(11));
t47 = cos(pkin(11));
t12 = t45 * t23 + t47 * t27;
t10 = qJD(3) * pkin(8) + t12;
t48 = cos(pkin(10));
t55 = -pkin(1) * t48 - pkin(2);
t32 = qJD(4) + (-pkin(3) * t52 + t55) * qJD(1);
t58 = qJD(1) * t52;
t59 = qJD(1) * t51;
t33 = t45 * t59 - t47 * t58;
t35 = (t45 * t52 + t47 * t51) * qJD(1);
t18 = t33 * pkin(4) - t35 * pkin(8) + t32;
t50 = sin(qJ(5));
t61 = cos(qJ(5));
t6 = t61 * t10 + t50 * t18;
t60 = cos(qJ(6));
t56 = qJD(1) * qJD(3);
t5 = -t50 * t10 + t61 * t18;
t11 = t47 * t23 - t45 * t27;
t31 = qJD(5) + t33;
t9 = -qJD(3) * pkin(4) - t11;
t49 = sin(qJ(6));
t43 = qJD(3) ^ 2 / 0.2e1;
t38 = t55 * qJD(1);
t29 = -t51 * t37 + t42;
t28 = qJD(6) + t31;
t26 = t50 * qJD(3) + t61 * t35;
t24 = -t61 * qJD(3) + t50 * t35;
t15 = -t49 * t24 + t60 * t26;
t13 = t60 * t24 + t49 * t26;
t7 = t24 * pkin(5) + t9;
t4 = -t24 * pkin(9) + t6;
t3 = t31 * pkin(5) - t26 * pkin(9) + t5;
t2 = t49 * t3 + t60 * t4;
t1 = t60 * t3 - t49 * t4;
t8 = [0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t48 * t62, -t46 * t62, 0, qJD(2) ^ 2 / 0.2e1 + (t46 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t53, t51 ^ 2 * t44, t51 * t53 * t52, t51 * t56, t52 ^ 2 * t44, t52 * t56, t43, t29 * qJD(3) - t38 * t58, -t30 * qJD(3) + t38 * t59 (-t29 * t51 + t30 * t52) * qJD(1), t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, qJD(3) * t35, t33 ^ 2 / 0.2e1, -qJD(3) * t33, t43, t11 * qJD(3) + t32 * t33, -t12 * qJD(3) + t32 * t35, -t11 * t35 - t12 * t33, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1, t26 ^ 2 / 0.2e1, -t26 * t24, t26 * t31, t24 ^ 2 / 0.2e1, -t24 * t31, t31 ^ 2 / 0.2e1, t9 * t24 + t5 * t31, t9 * t26 - t6 * t31, -t6 * t24 - t5 * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t28, t13 ^ 2 / 0.2e1, -t13 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t13, t7 * t15 - t2 * t28, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
