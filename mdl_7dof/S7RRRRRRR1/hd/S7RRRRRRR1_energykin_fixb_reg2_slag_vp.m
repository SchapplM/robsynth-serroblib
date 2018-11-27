% Calculate inertial parameters regressor of fixed base kinetic energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% T_reg [1x(7*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S7RRRRRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 20:34:12
% EndTime: 2018-11-26 20:34:12
% DurationCPUTime: 0.26s
% Computational Cost: add. (1024->57), mult. (2170->144), div. (0->0), fcn. (1795->12), ass. (0->51)
t48 = sin(qJ(3));
t52 = cos(qJ(3));
t49 = sin(qJ(2));
t60 = qJD(1) * t49;
t41 = -t48 * qJD(2) + t52 * t60;
t40 = -t52 * qJD(2) - t48 * t60;
t37 = t40 * pkin(2);
t66 = t37 ^ 2 / 0.2e1;
t54 = qJD(1) ^ 2;
t65 = t54 / 0.2e1;
t64 = cos(qJ(5));
t63 = cos(qJ(6));
t47 = sin(qJ(4));
t62 = t37 * t47;
t51 = cos(qJ(4));
t61 = t51 * t37;
t38 = -qJD(4) - t40;
t25 = -t38 * pkin(3) + t61;
t46 = sin(qJ(5));
t53 = cos(qJ(2));
t42 = t53 * qJD(1) + qJD(3);
t30 = t51 * t41 - t47 * t42;
t36 = t41 * pkin(2);
t55 = -t30 * pkin(3) - t36;
t16 = t64 * t25 + t46 * t55;
t45 = sin(qJ(6));
t10 = t63 * t16 + t45 * t62;
t14 = t46 * t25 - t64 * t55;
t58 = t14 ^ 2 / 0.2e1;
t57 = qJD(1) * qJD(2);
t23 = t64 * t30 - t46 * t38;
t28 = t47 * t41 + t51 * t42;
t27 = qJD(5) + t28;
t17 = t45 * t23 - t63 * t27;
t21 = t46 * t30 + t64 * t38;
t50 = cos(qJ(7));
t44 = sin(qJ(7));
t34 = t36 ^ 2 / 0.2e1;
t32 = t47 ^ 2 * t66;
t20 = qJD(6) + t21;
t19 = t63 * t23 + t45 * t27;
t12 = -qJD(7) + t17;
t9 = -t45 * t16 + t63 * t62;
t8 = t9 ^ 2 / 0.2e1;
t7 = t50 * t19 - t44 * t20;
t5 = t44 * t19 + t50 * t20;
t4 = -t20 * pkin(4) + t10;
t3 = t19 * pkin(4) + t14;
t2 = -t44 * t3 + t50 * t4;
t1 = -t50 * t3 - t44 * t4;
t6 = [0, 0, 0, 0, 0, t65, 0, 0, 0, 0, t49 ^ 2 * t65, t49 * t54 * t53, t49 * t57, t53 ^ 2 * t65, t53 * t57, qJD(2) ^ 2 / 0.2e1, 0, 0, 0, 0, t41 ^ 2 / 0.2e1, t41 * t40, t41 * t42, t40 ^ 2 / 0.2e1, t40 * t42, t42 ^ 2 / 0.2e1, -t36 * t42, -t37 * t42, t36 * t41 + t37 * t40, t66 + t34, t30 ^ 2 / 0.2e1, -t30 * t28, -t30 * t38, t28 ^ 2 / 0.2e1, t28 * t38, t38 ^ 2 / 0.2e1, -t36 * t28 + t38 * t62, -t36 * t30 + t38 * t61 (-t28 * t51 + t30 * t47) * t37, t51 ^ 2 * t66 + t32 + t34, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * t27, t21 ^ 2 / 0.2e1, -t21 * t27, t27 ^ 2 / 0.2e1, -t14 * t27 + t21 * t62, -t16 * t27 + t23 * t62, t14 * t23 - t16 * t21, t16 ^ 2 / 0.2e1 + t58 + t32, t19 ^ 2 / 0.2e1, -t19 * t17, t19 * t20, t17 ^ 2 / 0.2e1, -t17 * t20, t20 ^ 2 / 0.2e1, t14 * t17 + t9 * t20, -t10 * t20 + t14 * t19, -t10 * t17 - t9 * t19, t10 ^ 2 / 0.2e1 + t8 + t58, t7 ^ 2 / 0.2e1, -t7 * t5, -t7 * t12, t5 ^ 2 / 0.2e1, t5 * t12, t12 ^ 2 / 0.2e1, -t1 * t12 + t9 * t5, t2 * t12 + t9 * t7, -t1 * t7 - t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8;];
T_reg  = t6;
