% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:33
% EndTime: 2019-03-09 05:13:33
% DurationCPUTime: 0.20s
% Computational Cost: add. (690->64), mult. (1790->135), div. (0->0), fcn. (1323->8), ass. (0->52)
t50 = qJD(1) ^ 2;
t65 = t50 / 0.2e1;
t64 = pkin(4) + pkin(9);
t63 = cos(qJ(3));
t62 = cos(qJ(6));
t48 = sin(qJ(3));
t45 = cos(pkin(10));
t56 = qJD(1) * t45;
t44 = sin(pkin(10));
t57 = qJD(1) * t44;
t32 = t48 * t57 - t63 * t56;
t34 = (t63 * t44 + t45 * t48) * qJD(1);
t47 = sin(qJ(4));
t49 = cos(qJ(4));
t21 = t49 * t32 + t47 * t34;
t23 = -t47 * t32 + t49 * t34;
t61 = t23 * t21;
t43 = qJD(3) + qJD(4);
t60 = t23 * t43;
t59 = t43 * t21;
t58 = pkin(7) + qJ(2);
t35 = t58 * t57;
t36 = t58 * t56;
t25 = -t63 * t35 - t48 * t36;
t14 = qJD(3) * pkin(3) - t34 * pkin(8) + t25;
t26 = -t48 * t35 + t63 * t36;
t15 = -t32 * pkin(8) + t26;
t9 = t47 * t14 + t49 * t15;
t55 = t21 ^ 2 / 0.2e1;
t54 = t23 ^ 2 / 0.2e1;
t8 = t49 * t14 - t47 * t15;
t7 = -t43 * qJ(5) - t9;
t53 = qJD(5) - t8;
t37 = qJD(2) + (-pkin(2) * t45 - pkin(1)) * qJD(1);
t27 = t32 * pkin(3) + t37;
t52 = -t23 * qJ(5) + t27;
t46 = sin(qJ(6));
t42 = t45 ^ 2;
t41 = t44 ^ 2;
t40 = -qJD(1) * pkin(1) + qJD(2);
t39 = t43 ^ 2 / 0.2e1;
t20 = qJD(6) + t23;
t18 = t46 * t21 + t62 * t43;
t16 = -t62 * t21 + t46 * t43;
t10 = t21 * pkin(4) + t52;
t6 = -t43 * pkin(4) + t53;
t5 = t64 * t21 + t52;
t4 = -t21 * pkin(5) - t7;
t3 = t23 * pkin(5) - t64 * t43 + t53;
t2 = t46 * t3 + t62 * t5;
t1 = t62 * t3 - t46 * t5;
t11 = [0, 0, 0, 0, 0, t65, 0, 0, 0, 0, t41 * t65, t44 * t50 * t45, 0, t42 * t65, 0, 0, -t40 * t56, t40 * t57 (t41 + t42) * t50 * qJ(2), t40 ^ 2 / 0.2e1 + (t42 / 0.2e1 + t41 / 0.2e1) * qJ(2) ^ 2 * t50, t34 ^ 2 / 0.2e1, -t34 * t32, t34 * qJD(3), t32 ^ 2 / 0.2e1, -t32 * qJD(3), qJD(3) ^ 2 / 0.2e1, t25 * qJD(3) + t37 * t32, -t26 * qJD(3) + t37 * t34, -t25 * t34 - t26 * t32, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t54, -t61, t60, t55, -t59, t39, t27 * t21 + t8 * t43, t27 * t23 - t9 * t43, -t9 * t21 - t8 * t23, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t39, -t60, t59, t54, -t61, t55, t7 * t21 + t6 * t23, -t10 * t21 + t6 * t43, -t10 * t23 - t7 * t43, t10 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t20, t16 ^ 2 / 0.2e1, -t16 * t20, t20 ^ 2 / 0.2e1, t1 * t20 + t4 * t16, t4 * t18 - t2 * t20, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg  = t11;
