% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:42
% EndTime: 2019-12-31 17:53:42
% DurationCPUTime: 0.14s
% Computational Cost: add. (147->38), mult. (412->76), div. (0->0), fcn. (221->4), ass. (0->37)
t36 = qJD(1) ^ 2;
t46 = t36 / 0.2e1;
t31 = sin(pkin(7));
t43 = qJD(1) * t31;
t17 = qJ(2) * t43 + qJD(3);
t15 = -pkin(6) * t43 + t17;
t32 = cos(pkin(7));
t42 = qJD(1) * t32;
t16 = (-pkin(6) + qJ(2)) * t42;
t34 = sin(qJ(4));
t35 = cos(qJ(4));
t5 = t34 * t15 + t35 * t16;
t12 = (-t31 * t34 - t32 * t35) * qJD(1);
t14 = -t34 * t42 + t35 * t43;
t45 = t14 * t12;
t44 = qJ(2) * t36;
t41 = qJD(4) * t12;
t40 = t12 ^ 2 / 0.2e1;
t27 = -qJD(1) * pkin(1) + qJD(2);
t39 = t31 * t36 * t32;
t38 = qJ(2) ^ 2 * t46;
t10 = -pkin(2) * t42 - qJ(3) * t43 + t27;
t6 = pkin(3) * t42 - t10;
t4 = t35 * t15 - t34 * t16;
t30 = qJD(4) ^ 2 / 0.2e1;
t29 = t32 ^ 2;
t28 = t31 ^ 2;
t24 = t29 * t44;
t21 = t29 * t46;
t20 = t28 * t46;
t18 = t29 * t38;
t8 = t14 * qJD(4);
t7 = t14 ^ 2 / 0.2e1;
t3 = qJD(4) * qJ(5) + t5;
t2 = -qJD(4) * pkin(4) + qJD(5) - t4;
t1 = -t12 * pkin(4) - t14 * qJ(5) + t6;
t9 = [0, 0, 0, 0, 0, t46, 0, 0, 0, 0, t20, t39, 0, t21, 0, 0, -t27 * t42, t27 * t43, t28 * t44 + t24, t18 + t28 * t38 + t27 ^ 2 / 0.2e1, t20, 0, -t39, 0, 0, t21, -t10 * t42, t17 * t43 + t24, -t10 * t43, t18 + t10 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t7, t45, t8, t40, t41, t30, t4 * qJD(4) - t6 * t12, -t5 * qJD(4) + t6 * t14, t5 * t12 - t4 * t14, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1, t7, t8, -t45, t30, -t41, t40, -t2 * qJD(4) - t1 * t12, t3 * t12 + t2 * t14, t3 * qJD(4) - t1 * t14, t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t9;
