% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:44
% EndTime: 2019-12-31 20:55:45
% DurationCPUTime: 0.17s
% Computational Cost: add. (375->46), mult. (973->108), div. (0->0), fcn. (647->6), ass. (0->42)
t37 = qJD(1) ^ 2;
t51 = t37 / 0.2e1;
t50 = -pkin(7) - pkin(6);
t33 = sin(pkin(8));
t45 = cos(pkin(8));
t35 = sin(qJ(2));
t44 = qJD(1) * t35;
t25 = qJD(2) * pkin(2) + t50 * t44;
t36 = cos(qJ(2));
t43 = qJD(1) * t36;
t26 = t50 * t43;
t34 = sin(qJ(3));
t49 = cos(qJ(3));
t15 = t49 * t25 + t34 * t26;
t23 = (t34 * t36 + t49 * t35) * qJD(1);
t30 = qJD(2) + qJD(3);
t7 = t30 * pkin(3) - t23 * qJ(4) + t15;
t16 = t34 * t25 - t49 * t26;
t21 = t34 * t44 - t49 * t43;
t9 = -t21 * qJ(4) + t16;
t4 = t33 * t7 + t45 * t9;
t12 = t45 * t21 + t33 * t23;
t14 = -t33 * t21 + t45 * t23;
t48 = t14 * t12;
t47 = t30 * t12;
t46 = t36 * t37;
t42 = t12 ^ 2 / 0.2e1;
t41 = qJD(1) * qJD(2);
t40 = t35 * t41;
t39 = t36 * t41;
t27 = (-pkin(2) * t36 - pkin(1)) * qJD(1);
t3 = -t33 * t9 + t45 * t7;
t17 = t21 * pkin(3) + qJD(4) + t27;
t32 = t36 ^ 2;
t31 = t35 ^ 2;
t29 = t30 ^ 2 / 0.2e1;
t11 = t14 ^ 2 / 0.2e1;
t10 = t14 * t30;
t5 = t12 * pkin(4) - t14 * qJ(5) + t17;
t2 = t30 * qJ(5) + t4;
t1 = -t30 * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t51, 0, 0, 0, 0, t31 * t51, t35 * t46, t40, t32 * t51, t39, qJD(2) ^ 2 / 0.2e1, pkin(1) * t46 - pkin(6) * t40, -t37 * pkin(1) * t35 - pkin(6) * t39, (t31 + t32) * t37 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t32 / 0.2e1 + t31 / 0.2e1) * pkin(6) ^ 2) * t37, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * t30, t21 ^ 2 / 0.2e1, -t21 * t30, t29, t15 * t30 + t27 * t21, -t16 * t30 + t27 * t23, -t15 * t23 - t16 * t21, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t11, -t48, t10, t42, -t47, t29, t17 * t12 + t3 * t30, t17 * t14 - t4 * t30, -t4 * t12 - t3 * t14, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t11, t10, t48, t29, t47, t42, -t1 * t30 + t5 * t12, t1 * t14 - t2 * t12, -t5 * t14 + t2 * t30, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
