% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:50
% EndTime: 2019-12-31 17:30:50
% DurationCPUTime: 0.16s
% Computational Cost: add. (283->40), mult. (766->104), div. (0->0), fcn. (539->8), ass. (0->36)
t47 = cos(qJ(3));
t46 = cos(qJ(4));
t31 = sin(pkin(4));
t37 = qJD(1) ^ 2;
t45 = t31 ^ 2 * t37;
t35 = sin(qJ(2));
t36 = cos(qJ(2));
t44 = qJD(1) * t31;
t39 = t36 * t44;
t43 = cos(pkin(4)) * qJD(1);
t41 = pkin(1) * t43;
t21 = pkin(6) * t39 + t35 * t41;
t29 = qJD(2) + t43;
t14 = t29 * pkin(7) + t21;
t15 = (-pkin(2) * t36 - pkin(7) * t35 - pkin(1)) * t44;
t34 = sin(qJ(3));
t7 = t47 * t14 + t34 * t15;
t42 = t36 * t45;
t40 = t35 * t44;
t38 = t45 / 0.2e1;
t20 = -pkin(6) * t40 + t36 * t41;
t17 = -t47 * t29 + t34 * t40;
t6 = -t34 * t14 + t47 * t15;
t13 = -t29 * pkin(2) - t20;
t33 = sin(qJ(4));
t24 = -qJD(3) + t39;
t19 = t34 * t29 + t47 * t40;
t16 = qJD(4) + t17;
t10 = t46 * t19 - t33 * t24;
t8 = t33 * t19 + t46 * t24;
t5 = -t24 * pkin(8) + t7;
t4 = t24 * pkin(3) - t6;
t3 = t17 * pkin(3) - t19 * pkin(8) + t13;
t2 = t33 * t3 + t46 * t5;
t1 = t46 * t3 - t33 * t5;
t9 = [0, 0, 0, 0, 0, t37 / 0.2e1, 0, 0, 0, 0, t35 ^ 2 * t38, t35 * t42, t29 * t40, t36 ^ 2 * t38, t29 * t39, t29 ^ 2 / 0.2e1, pkin(1) * t42 + t20 * t29, -pkin(1) * t35 * t45 - t21 * t29, (-t20 * t35 + t21 * t36) * t44, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t38, t19 ^ 2 / 0.2e1, -t19 * t17, -t19 * t24, t17 ^ 2 / 0.2e1, t17 * t24, t24 ^ 2 / 0.2e1, t13 * t17 - t6 * t24, t13 * t19 + t7 * t24, -t7 * t17 - t6 * t19, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t16, t8 ^ 2 / 0.2e1, -t8 * t16, t16 ^ 2 / 0.2e1, t1 * t16 + t4 * t8, t4 * t10 - t2 * t16, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg = t9;
