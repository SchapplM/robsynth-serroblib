% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:44
% EndTime: 2019-03-09 03:19:44
% DurationCPUTime: 0.18s
% Computational Cost: add. (455->60), mult. (1158->115), div. (0->0), fcn. (778->6), ass. (0->51)
t44 = sin(pkin(9));
t45 = cos(pkin(9));
t47 = sin(qJ(3));
t48 = cos(qJ(3));
t33 = (t44 * t48 + t45 * t47) * qJD(1);
t49 = qJD(1) ^ 2;
t64 = t49 / 0.2e1;
t63 = pkin(3) + pkin(8);
t59 = qJD(1) * t44;
t60 = pkin(7) + qJ(2);
t35 = t60 * t59;
t58 = qJD(1) * t45;
t36 = t60 * t58;
t18 = -t48 * t35 - t47 * t36;
t53 = qJD(4) - t18;
t10 = t33 * pkin(4) - t63 * qJD(3) + t53;
t46 = sin(qJ(5));
t62 = cos(qJ(5));
t31 = t47 * t59 - t48 * t58;
t37 = qJD(2) + (-pkin(2) * t45 - pkin(1)) * qJD(1);
t51 = -t33 * qJ(4) + t37;
t7 = t63 * t31 + t51;
t4 = t46 * t10 + t62 * t7;
t61 = t33 * t31;
t19 = -t47 * t35 + t48 * t36;
t57 = qJD(3) * t31;
t56 = t33 * qJD(3);
t55 = t31 ^ 2 / 0.2e1;
t54 = t33 ^ 2 / 0.2e1;
t17 = -qJD(3) * qJ(4) - t19;
t3 = t62 * t10 - t46 * t7;
t11 = -t31 * pkin(4) - t17;
t42 = qJD(3) ^ 2 / 0.2e1;
t41 = t45 ^ 2;
t40 = t44 ^ 2;
t39 = -qJD(1) * pkin(1) + qJD(2);
t27 = qJD(5) + t33;
t25 = t27 ^ 2 / 0.2e1;
t24 = t62 * qJD(3) + t46 * t31;
t22 = t46 * qJD(3) - t62 * t31;
t21 = t24 ^ 2 / 0.2e1;
t20 = t22 ^ 2 / 0.2e1;
t16 = -qJD(3) * pkin(3) + t53;
t15 = t24 * t27;
t14 = t22 * t27;
t13 = t31 * pkin(3) + t51;
t12 = t24 * t22;
t5 = t22 * pkin(5) + qJD(6) + t11;
t2 = -t22 * qJ(6) + t4;
t1 = t27 * pkin(5) - t24 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t64, 0, 0, 0, 0, t40 * t64, t44 * t49 * t45, 0, t41 * t64, 0, 0, -t39 * t58, t39 * t59 (t40 + t41) * t49 * qJ(2), t39 ^ 2 / 0.2e1 + (t41 / 0.2e1 + t40 / 0.2e1) * qJ(2) ^ 2 * t49, t54, -t61, t56, t55, -t57, t42, t18 * qJD(3) + t37 * t31, -t19 * qJD(3) + t37 * t33, -t18 * t33 - t19 * t31, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t42, -t56, t57, t54, -t61, t55, t16 * t33 + t17 * t31, t16 * qJD(3) - t13 * t31, -t17 * qJD(3) - t13 * t33, t13 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t21, -t12, t15, t20, -t14, t25, t11 * t22 + t3 * t27, t11 * t24 - t4 * t27, -t4 * t22 - t3 * t24, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1, t21, -t12, t15, t20, -t14, t25, t1 * t27 + t5 * t22, -t2 * t27 + t5 * t24, -t1 * t24 - t2 * t22, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;
