% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:23
% EndTime: 2019-12-31 21:09:23
% DurationCPUTime: 0.16s
% Computational Cost: add. (227->46), mult. (544->94), div. (0->0), fcn. (296->4), ass. (0->38)
t33 = qJD(1) ^ 2;
t46 = t33 / 0.2e1;
t45 = cos(qJ(3));
t30 = sin(qJ(3));
t31 = sin(qJ(2));
t41 = qJD(1) * t31;
t15 = -t45 * qJD(2) + t30 * t41;
t17 = t30 * qJD(2) + t45 * t41;
t44 = t15 * t17;
t32 = cos(qJ(2));
t40 = t32 * qJD(1);
t24 = -qJD(3) + t40;
t9 = t17 * t24;
t10 = t24 * t15;
t43 = t32 * t33;
t42 = pkin(3) + qJ(5);
t12 = (-pkin(2) * t32 - pkin(7) * t31 - pkin(1)) * qJD(1);
t22 = pkin(6) * t40 + qJD(2) * pkin(7);
t8 = t30 * t12 + t45 * t22;
t13 = t15 ^ 2 / 0.2e1;
t14 = t17 ^ 2 / 0.2e1;
t39 = qJD(1) * qJD(2);
t38 = t31 * t39;
t37 = t32 * t39;
t21 = -qJD(2) * pkin(2) + pkin(6) * t41;
t5 = t24 * qJ(4) - t8;
t7 = t45 * t12 - t30 * t22;
t36 = qJD(4) - t7;
t35 = -t17 * qJ(4) + t21;
t29 = t32 ^ 2;
t28 = t31 ^ 2;
t23 = t24 ^ 2 / 0.2e1;
t6 = t15 * pkin(3) + t35;
t4 = t24 * pkin(3) + t36;
t3 = t42 * t15 + t35;
t2 = -t15 * pkin(4) + qJD(5) - t5;
t1 = t17 * pkin(4) + t42 * t24 + t36;
t11 = [0, 0, 0, 0, 0, t46, 0, 0, 0, 0, t28 * t46, t31 * t43, t38, t29 * t46, t37, qJD(2) ^ 2 / 0.2e1, pkin(1) * t43 - pkin(6) * t38, -t33 * pkin(1) * t31 - pkin(6) * t37, (t28 + t29) * t33 * pkin(6), (pkin(1) ^ 2 / 0.2e1 + (t29 / 0.2e1 + t28 / 0.2e1) * pkin(6) ^ 2) * t33, t14, -t44, -t9, t13, t10, t23, t21 * t15 - t7 * t24, t21 * t17 + t8 * t24, -t8 * t15 - t7 * t17, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t23, t9, -t10, t14, -t44, t13, t5 * t15 + t4 * t17, -t6 * t15 - t4 * t24, -t6 * t17 + t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t23, -t10, -t9, t13, t44, t14, t1 * t17 - t2 * t15, -t3 * t17 - t2 * t24, t1 * t24 + t3 * t15, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t11;
