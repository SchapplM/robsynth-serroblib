% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:39
% EndTime: 2019-03-08 22:52:39
% DurationCPUTime: 0.16s
% Computational Cost: add. (317->52), mult. (737->115), div. (0->0), fcn. (481->8), ass. (0->46)
t43 = qJD(2) ^ 2;
t58 = t43 / 0.2e1;
t57 = pkin(4) + pkin(5);
t56 = cos(qJ(4));
t38 = sin(qJ(4));
t39 = sin(qJ(3));
t52 = qJD(2) * t39;
t23 = -t56 * qJD(3) + t38 * t52;
t41 = cos(qJ(3));
t51 = t41 * qJD(2);
t32 = -qJD(4) + t51;
t55 = t23 * t32;
t25 = t38 * qJD(3) + t56 * t52;
t11 = t25 * t23;
t18 = t25 * t32;
t40 = sin(qJ(2));
t36 = sin(pkin(6));
t54 = qJD(1) * t36;
t26 = qJD(2) * pkin(8) + t40 * t54;
t37 = cos(pkin(6));
t53 = qJD(1) * t37;
t16 = t41 * t26 + t39 * t53;
t13 = qJD(3) * pkin(9) + t16;
t42 = cos(qJ(2));
t49 = t42 * t54;
t17 = -t49 + (-pkin(3) * t41 - pkin(9) * t39 - pkin(2)) * qJD(2);
t7 = t56 * t13 + t38 * t17;
t15 = -t39 * t26 + t41 * t53;
t19 = t23 ^ 2 / 0.2e1;
t28 = t32 ^ 2 / 0.2e1;
t50 = qJD(2) * qJD(3);
t5 = -t32 * qJ(5) + t7;
t48 = qJD(2) * t54;
t47 = qJD(3) * pkin(3) + t15;
t6 = -t38 * t13 + t56 * t17;
t46 = qJD(5) - t6;
t45 = t25 * qJ(5) + t47;
t44 = qJD(1) ^ 2;
t27 = -qJD(2) * pkin(2) - t49;
t20 = t25 ^ 2 / 0.2e1;
t8 = t23 * pkin(4) - t45;
t4 = t32 * pkin(4) + t46;
t3 = -t57 * t23 + qJD(6) + t45;
t2 = t23 * qJ(6) + t5;
t1 = -t25 * qJ(6) + t57 * t32 + t46;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t44 / 0.2e1, 0, 0, 0, 0, 0, t58, t42 * t48, -t40 * t48, 0 (t37 ^ 2 / 0.2e1 + (t40 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1) * t36 ^ 2) * t44, t39 ^ 2 * t58, t39 * t43 * t41, t39 * t50, t41 ^ 2 * t58, t41 * t50, qJD(3) ^ 2 / 0.2e1, t15 * qJD(3) - t27 * t51, -t16 * qJD(3) + t27 * t52 (-t15 * t39 + t16 * t41) * qJD(2), t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t20, -t11, -t18, t19, t55, t28, -t23 * t47 - t6 * t32, -t25 * t47 + t7 * t32, -t7 * t23 - t6 * t25, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1, t20, -t18, t11, t28, -t55, t19, t8 * t23 + t4 * t32, -t5 * t23 + t4 * t25, -t8 * t25 - t5 * t32, t5 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t20, t11, t18, t19, t55, t28, t1 * t32 - t3 * t23, -t2 * t32 + t3 * t25, -t1 * t25 + t2 * t23, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t9;
