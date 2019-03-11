% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:36
% EndTime: 2019-03-09 02:51:36
% DurationCPUTime: 0.22s
% Computational Cost: add. (612->64), mult. (1529->133), div. (0->0), fcn. (1068->8), ass. (0->51)
t46 = sin(pkin(9));
t47 = cos(pkin(9));
t49 = sin(qJ(3));
t50 = cos(qJ(3));
t34 = (t46 * t50 + t47 * t49) * qJD(1);
t51 = qJD(1) ^ 2;
t66 = t51 / 0.2e1;
t60 = qJD(1) * t46;
t62 = pkin(7) + qJ(2);
t36 = t62 * t60;
t59 = qJD(1) * t47;
t37 = t62 * t59;
t20 = -t50 * t36 - t49 * t37;
t55 = qJD(4) - t20;
t63 = pkin(3) + qJ(5);
t12 = t34 * pkin(4) - t63 * qJD(3) + t55;
t45 = sin(pkin(10));
t61 = cos(pkin(10));
t32 = t49 * t60 - t50 * t59;
t38 = qJD(2) + (-pkin(2) * t47 - pkin(1)) * qJD(1);
t53 = -t34 * qJ(4) + t38;
t9 = t63 * t32 + t53;
t6 = t45 * t12 + t61 * t9;
t65 = cos(qJ(6));
t64 = t34 * t32;
t21 = -t49 * t36 + t50 * t37;
t58 = qJD(3) * t32;
t57 = t34 * qJD(3);
t56 = t32 ^ 2 / 0.2e1;
t27 = t34 ^ 2 / 0.2e1;
t19 = -qJD(3) * qJ(4) - t21;
t5 = t61 * t12 - t45 * t9;
t16 = -t32 * pkin(4) + qJD(5) - t19;
t48 = sin(qJ(6));
t43 = qJD(3) ^ 2 / 0.2e1;
t42 = t47 ^ 2;
t41 = t46 ^ 2;
t40 = -qJD(1) * pkin(1) + qJD(2);
t28 = qJD(6) + t34;
t25 = t61 * qJD(3) + t45 * t32;
t23 = t45 * qJD(3) - t61 * t32;
t18 = -qJD(3) * pkin(3) + t55;
t17 = t32 * pkin(3) + t53;
t15 = -t48 * t23 + t65 * t25;
t13 = t65 * t23 + t48 * t25;
t7 = t23 * pkin(5) + t16;
t4 = -t23 * pkin(8) + t6;
t3 = t34 * pkin(5) - t25 * pkin(8) + t5;
t2 = t48 * t3 + t65 * t4;
t1 = t65 * t3 - t48 * t4;
t8 = [0, 0, 0, 0, 0, t66, 0, 0, 0, 0, t41 * t66, t46 * t51 * t47, 0, t42 * t66, 0, 0, -t40 * t59, t40 * t60 (t41 + t42) * t51 * qJ(2), t40 ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * qJ(2) ^ 2 * t51, t27, -t64, t57, t56, -t58, t43, t20 * qJD(3) + t38 * t32, -t21 * qJD(3) + t38 * t34, -t20 * t34 - t21 * t32, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t38 ^ 2 / 0.2e1, t43, -t57, t58, t27, -t64, t56, t18 * t34 + t19 * t32, t18 * qJD(3) - t17 * t32, -t19 * qJD(3) - t17 * t34, t17 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t25 ^ 2 / 0.2e1, -t25 * t23, t25 * t34, t23 ^ 2 / 0.2e1, -t23 * t34, t27, t16 * t23 + t5 * t34, t16 * t25 - t6 * t34, -t6 * t23 - t5 * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t28, t13 ^ 2 / 0.2e1, -t13 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t7 * t13, t7 * t15 - t2 * t28, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
