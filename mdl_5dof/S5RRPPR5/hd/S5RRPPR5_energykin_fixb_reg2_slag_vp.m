% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:56
% EndTime: 2019-12-31 19:29:57
% DurationCPUTime: 0.20s
% Computational Cost: add. (282->47), mult. (760->108), div. (0->0), fcn. (467->6), ass. (0->45)
t33 = sin(pkin(8));
t35 = sin(qJ(2));
t36 = cos(qJ(2));
t48 = cos(pkin(8));
t20 = (t33 * t36 + t48 * t35) * qJD(1);
t25 = -(pkin(2) * t36 + pkin(1)) * qJD(1) + qJD(3);
t56 = -t20 * qJ(4) + t25;
t37 = qJD(1) ^ 2;
t55 = t37 / 0.2e1;
t54 = -pkin(3) - pkin(4);
t53 = cos(qJ(5));
t46 = qJD(1) * t36;
t47 = qJD(1) * t35;
t18 = t33 * t47 - t48 * t46;
t52 = t20 * t18;
t51 = t36 * t37;
t50 = pkin(6) + qJ(3);
t23 = qJD(2) * pkin(2) - t50 * t47;
t24 = t50 * t46;
t13 = t33 * t23 + t48 * t24;
t45 = qJD(2) * t18;
t44 = t18 ^ 2 / 0.2e1;
t43 = qJD(1) * qJD(2);
t11 = qJD(2) * qJ(4) + t13;
t41 = t35 * t43;
t40 = t36 * t43;
t12 = t48 * t23 - t33 * t24;
t39 = qJD(4) - t12;
t34 = sin(qJ(5));
t32 = t36 ^ 2;
t31 = t35 ^ 2;
t29 = qJD(2) ^ 2 / 0.2e1;
t27 = qJD(2) - qJD(5);
t17 = t20 * qJD(2);
t16 = t20 ^ 2 / 0.2e1;
t10 = t34 * t18 + t53 * t20;
t8 = -t53 * t18 + t34 * t20;
t7 = -qJD(2) * pkin(3) + t39;
t6 = t18 * pkin(3) + t56;
t5 = t18 * pkin(7) + t11;
t4 = -t20 * pkin(7) + t54 * qJD(2) + t39;
t3 = t54 * t18 - t56;
t2 = t34 * t4 + t53 * t5;
t1 = -t34 * t5 + t53 * t4;
t9 = [0, 0, 0, 0, 0, t55, 0, 0, 0, 0, t31 * t55, t35 * t51, t41, t32 * t55, t40, t29, pkin(1) * t51 - pkin(6) * t41, -t37 * pkin(1) * t35 - pkin(6) * t40, (t31 + t32) * t37 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t32 / 0.2e1 + t31 / 0.2e1) * pkin(6) ^ 2) * t37, t16, -t52, t17, t44, -t45, t29, t12 * qJD(2) + t25 * t18, -t13 * qJD(2) + t25 * t20, -t12 * t20 - t13 * t18, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t16, t17, t52, t29, t45, t44, -t7 * qJD(2) + t6 * t18, -t11 * t18 + t7 * t20, t11 * qJD(2) - t6 * t20, t11 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, -t10 * t27, t8 ^ 2 / 0.2e1, t8 * t27, t27 ^ 2 / 0.2e1, -t1 * t27 + t3 * t8, t3 * t10 + t2 * t27, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t9;
