% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPR7
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:21
% EndTime: 2019-12-31 19:36:21
% DurationCPUTime: 0.19s
% Computational Cost: add. (280->49), mult. (743->108), div. (0->0), fcn. (456->6), ass. (0->45)
t30 = sin(pkin(8));
t31 = cos(pkin(8));
t33 = sin(qJ(2));
t34 = cos(qJ(2));
t19 = (t30 * t34 + t31 * t33) * qJD(1);
t35 = qJD(1) ^ 2;
t54 = t35 / 0.2e1;
t53 = pkin(3) + pkin(7);
t52 = cos(qJ(5));
t47 = qJD(1) * t34;
t48 = qJD(1) * t33;
t17 = t30 * t48 - t31 * t47;
t51 = t19 * t17;
t50 = t34 * t35;
t49 = pkin(6) + qJ(3);
t23 = qJD(2) * pkin(2) - t49 * t48;
t24 = t49 * t47;
t10 = t30 * t23 + t31 * t24;
t46 = qJD(2) * t17;
t45 = t19 * qJD(2);
t44 = t17 ^ 2 / 0.2e1;
t43 = t19 ^ 2 / 0.2e1;
t42 = qJD(1) * qJD(2);
t41 = t33 * t42;
t40 = t34 * t42;
t9 = t31 * t23 - t30 * t24;
t39 = qJD(4) - t9;
t8 = -qJD(2) * qJ(4) - t10;
t25 = qJD(3) + (-pkin(2) * t34 - pkin(1)) * qJD(1);
t37 = -t19 * qJ(4) + t25;
t32 = sin(qJ(5));
t29 = t34 ^ 2;
t28 = t33 ^ 2;
t27 = qJD(2) ^ 2 / 0.2e1;
t16 = qJD(5) + t19;
t13 = t52 * qJD(2) + t32 * t17;
t11 = t32 * qJD(2) - t52 * t17;
t7 = -qJD(2) * pkin(3) + t39;
t6 = t17 * pkin(3) + t37;
t5 = -t17 * pkin(4) - t8;
t4 = t19 * pkin(4) - t53 * qJD(2) + t39;
t3 = t53 * t17 + t37;
t2 = t52 * t3 + t32 * t4;
t1 = -t32 * t3 + t52 * t4;
t12 = [0, 0, 0, 0, 0, t54, 0, 0, 0, 0, t28 * t54, t33 * t50, t41, t29 * t54, t40, t27, pkin(1) * t50 - pkin(6) * t41, -t35 * pkin(1) * t33 - pkin(6) * t40, (t28 + t29) * t35 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t29 / 0.2e1 + t28 / 0.2e1) * pkin(6) ^ 2) * t35, t43, -t51, t45, t44, -t46, t27, t9 * qJD(2) + t25 * t17, -t10 * qJD(2) + t25 * t19, -t10 * t17 - t9 * t19, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t27, -t45, t46, t43, -t51, t44, t8 * t17 + t7 * t19, t7 * qJD(2) - t6 * t17, -t8 * qJD(2) - t6 * t19, t6 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t16, t11 ^ 2 / 0.2e1, -t11 * t16, t16 ^ 2 / 0.2e1, t1 * t16 + t5 * t11, t5 * t13 - t2 * t16, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t12;
