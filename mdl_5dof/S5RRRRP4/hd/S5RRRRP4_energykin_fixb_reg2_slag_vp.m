% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:12
% EndTime: 2019-12-31 21:51:12
% DurationCPUTime: 0.13s
% Computational Cost: add. (300->36), mult. (458->93), div. (0->0), fcn. (239->6), ass. (0->41)
t27 = sin(qJ(3));
t24 = t27 ^ 2;
t46 = t24 / 0.2e1;
t29 = cos(qJ(3));
t25 = t29 ^ 2;
t45 = t25 / 0.2e1;
t23 = qJD(1) + qJD(2);
t28 = sin(qJ(2));
t39 = pkin(1) * qJD(1);
t35 = t28 * t39;
t17 = t23 * pkin(7) + t35;
t33 = pkin(8) * t23 + t17;
t10 = t33 * t29;
t26 = sin(qJ(4));
t44 = cos(qJ(4));
t8 = qJD(3) * pkin(3) - t33 * t27;
t5 = t44 * t10 + t26 * t8;
t40 = t23 * t29;
t41 = t23 * t27;
t12 = t26 * t41 - t44 * t40;
t14 = (t26 * t29 + t44 * t27) * t23;
t43 = t14 * t12;
t22 = qJD(3) + qJD(4);
t42 = t22 * t12;
t38 = qJD(3) * t27;
t37 = qJD(3) * t29;
t36 = t12 ^ 2 / 0.2e1;
t30 = cos(qJ(2));
t34 = t30 * t39;
t4 = -t26 * t10 + t44 * t8;
t15 = -t34 + (-pkin(3) * t29 - pkin(2)) * t23;
t31 = qJD(1) ^ 2;
t21 = t23 ^ 2;
t20 = t22 ^ 2 / 0.2e1;
t18 = -t23 * pkin(2) - t34;
t11 = t14 ^ 2 / 0.2e1;
t9 = t14 * t22;
t3 = t12 * pkin(4) - t14 * qJ(5) + t15;
t2 = t22 * qJ(5) + t5;
t1 = -t22 * pkin(4) + qJD(5) - t4;
t6 = [0, 0, 0, 0, 0, t31 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 / 0.2e1, t23 * t34, -t23 * t35, 0, (t28 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t31, t21 * t46, t27 * t21 * t29, t23 * t38, t21 * t45, t23 * t37, qJD(3) ^ 2 / 0.2e1, -t17 * t38 - t18 * t40, -t17 * t37 + t18 * t41, (t24 + t25) * t23 * t17, t18 ^ 2 / 0.2e1 + (t45 + t46) * t17 ^ 2, t11, -t43, t9, t36, -t42, t20, t15 * t12 + t4 * t22, t15 * t14 - t5 * t22, -t5 * t12 - t4 * t14, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t11, t9, t43, t20, t42, t36, -t1 * t22 + t3 * t12, t1 * t14 - t2 * t12, -t3 * t14 + t2 * t22, t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t6;
