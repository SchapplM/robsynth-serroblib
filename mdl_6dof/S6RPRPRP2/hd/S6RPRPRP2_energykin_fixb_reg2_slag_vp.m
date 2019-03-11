% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRP2
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
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:12
% EndTime: 2019-03-09 03:06:12
% DurationCPUTime: 0.16s
% Computational Cost: add. (462->55), mult. (1076->127), div. (0->0), fcn. (685->8), ass. (0->47)
t46 = qJD(1) ^ 2;
t39 = t46 / 0.2e1;
t58 = pkin(1) * t46;
t45 = cos(qJ(3));
t42 = cos(pkin(9));
t48 = -pkin(1) * t42 - pkin(2);
t27 = qJD(4) + (-pkin(3) * t45 + t48) * qJD(1);
t40 = sin(pkin(10));
t52 = qJD(1) * t45;
t44 = sin(qJ(3));
t53 = qJD(1) * t44;
t54 = cos(pkin(10));
t28 = t40 * t53 - t54 * t52;
t30 = (t40 * t45 + t54 * t44) * qJD(1);
t12 = t28 * pkin(4) - t30 * pkin(8) + t27;
t43 = sin(qJ(5));
t57 = cos(qJ(5));
t41 = sin(pkin(9));
t32 = (pkin(1) * t41 + pkin(7)) * qJD(1);
t37 = t45 * qJD(2);
t50 = qJ(4) * qJD(1);
t17 = qJD(3) * pkin(3) + t37 + (-t32 - t50) * t44;
t25 = t44 * qJD(2) + t45 * t32;
t22 = t45 * t50 + t25;
t10 = t40 * t17 + t54 * t22;
t8 = qJD(3) * pkin(8) + t10;
t4 = t43 * t12 + t57 * t8;
t19 = -t57 * qJD(3) + t43 * t30;
t21 = t43 * qJD(3) + t57 * t30;
t56 = t21 * t19;
t26 = qJD(5) + t28;
t55 = t26 * t19;
t51 = t19 ^ 2 / 0.2e1;
t49 = qJD(1) * qJD(3);
t9 = t54 * t17 - t40 * t22;
t3 = t57 * t12 - t43 * t8;
t7 = -qJD(3) * pkin(4) - t9;
t38 = qJD(3) ^ 2 / 0.2e1;
t33 = t48 * qJD(1);
t24 = -t44 * t32 + t37;
t23 = t26 ^ 2 / 0.2e1;
t18 = t21 ^ 2 / 0.2e1;
t13 = t21 * t26;
t5 = t19 * pkin(5) - t21 * qJ(6) + t7;
t2 = t26 * qJ(6) + t4;
t1 = -t26 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t42 * t58, -t41 * t58, 0, qJD(2) ^ 2 / 0.2e1 + (t41 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t46, t44 ^ 2 * t39, t44 * t46 * t45, t44 * t49, t45 ^ 2 * t39, t45 * t49, t38, t24 * qJD(3) - t33 * t52, -t25 * qJD(3) + t33 * t53 (-t24 * t44 + t25 * t45) * qJD(1), t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t28, t30 * qJD(3), t28 ^ 2 / 0.2e1, -t28 * qJD(3), t38, t9 * qJD(3) + t27 * t28, -t10 * qJD(3) + t27 * t30, -t10 * t28 - t9 * t30, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t18, -t56, t13, t51, -t55, t23, t7 * t19 + t3 * t26, t7 * t21 - t4 * t26, -t4 * t19 - t3 * t21, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t18, t13, t56, t23, t55, t51, -t1 * t26 + t5 * t19, t1 * t21 - t2 * t19, t2 * t26 - t5 * t21, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
