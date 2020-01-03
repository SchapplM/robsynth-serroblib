% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPP6
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:01:58
% EndTime: 2019-12-31 21:01:58
% DurationCPUTime: 0.17s
% Computational Cost: add. (375->49), mult. (884->109), div. (0->0), fcn. (558->6), ass. (0->41)
t39 = qJD(1) ^ 2;
t52 = t39 / 0.2e1;
t35 = sin(pkin(8));
t47 = cos(pkin(8));
t37 = sin(qJ(2));
t38 = cos(qJ(2));
t21 = (-pkin(2) * t38 - pkin(7) * t37 - pkin(1)) * qJD(1);
t45 = t38 * qJD(1);
t27 = pkin(6) * t45 + qJD(2) * pkin(7);
t36 = sin(qJ(3));
t51 = cos(qJ(3));
t15 = t51 * t21 - t36 * t27;
t46 = qJD(1) * t37;
t24 = t36 * qJD(2) + t51 * t46;
t29 = -qJD(3) + t45;
t7 = -t29 * pkin(3) - t24 * qJ(4) + t15;
t16 = t36 * t21 + t51 * t27;
t22 = -t51 * qJD(2) + t36 * t46;
t9 = -t22 * qJ(4) + t16;
t4 = t35 * t7 + t47 * t9;
t12 = t47 * t22 + t35 * t24;
t14 = -t35 * t22 + t47 * t24;
t50 = t14 * t12;
t49 = t29 * t12;
t48 = t38 * t39;
t44 = t12 ^ 2 / 0.2e1;
t43 = qJD(1) * qJD(2);
t42 = t37 * t43;
t41 = t38 * t43;
t26 = -qJD(2) * pkin(2) + pkin(6) * t46;
t3 = -t35 * t9 + t47 * t7;
t17 = t22 * pkin(3) + qJD(4) + t26;
t34 = t38 ^ 2;
t33 = t37 ^ 2;
t28 = t29 ^ 2 / 0.2e1;
t11 = t14 ^ 2 / 0.2e1;
t10 = t14 * t29;
t5 = t12 * pkin(4) - t14 * qJ(5) + t17;
t2 = -t29 * qJ(5) + t4;
t1 = t29 * pkin(4) + qJD(5) - t3;
t6 = [0, 0, 0, 0, 0, t52, 0, 0, 0, 0, t33 * t52, t37 * t48, t42, t34 * t52, t41, qJD(2) ^ 2 / 0.2e1, pkin(1) * t48 - pkin(6) * t42, -t39 * pkin(1) * t37 - pkin(6) * t41, (t33 + t34) * t39 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t34 / 0.2e1 + t33 / 0.2e1) * pkin(6) ^ 2) * t39, t24 ^ 2 / 0.2e1, -t24 * t22, -t24 * t29, t22 ^ 2 / 0.2e1, t22 * t29, t28, -t15 * t29 + t26 * t22, t16 * t29 + t26 * t24, -t15 * t24 - t16 * t22, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t11, -t50, -t10, t44, t49, t28, t17 * t12 - t3 * t29, t17 * t14 + t4 * t29, -t4 * t12 - t3 * t14, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t11, -t10, t50, t28, -t49, t44, t1 * t29 + t5 * t12, t1 * t14 - t2 * t12, -t5 * t14 - t2 * t29, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
