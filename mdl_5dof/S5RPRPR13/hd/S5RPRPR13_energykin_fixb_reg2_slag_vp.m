% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR13_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:51
% EndTime: 2019-12-31 18:32:52
% DurationCPUTime: 0.18s
% Computational Cost: add. (261->48), mult. (709->102), div. (0->0), fcn. (452->6), ass. (0->42)
t31 = sin(pkin(8));
t32 = cos(pkin(8));
t34 = sin(qJ(3));
t35 = cos(qJ(3));
t21 = (t31 * t35 + t32 * t34) * qJD(1);
t36 = qJD(1) ^ 2;
t51 = t36 / 0.2e1;
t50 = pkin(3) + pkin(7);
t49 = cos(qJ(5));
t45 = qJD(1) * t32;
t46 = qJD(1) * t31;
t19 = t34 * t46 - t35 * t45;
t48 = t21 * t19;
t47 = pkin(6) + qJ(2);
t23 = t47 * t46;
t24 = t47 * t45;
t10 = -t34 * t23 + t35 * t24;
t44 = qJD(3) * t19;
t43 = t21 * qJD(3);
t42 = t19 ^ 2 / 0.2e1;
t41 = t21 ^ 2 / 0.2e1;
t9 = -t35 * t23 - t34 * t24;
t40 = qJD(4) - t9;
t8 = -qJD(3) * qJ(4) - t10;
t25 = qJD(2) + (-pkin(2) * t32 - pkin(1)) * qJD(1);
t38 = -t21 * qJ(4) + t25;
t33 = sin(qJ(5));
t30 = qJD(3) ^ 2 / 0.2e1;
t29 = t32 ^ 2;
t28 = t31 ^ 2;
t27 = -qJD(1) * pkin(1) + qJD(2);
t15 = qJD(5) + t21;
t13 = t49 * qJD(3) + t33 * t19;
t11 = t33 * qJD(3) - t49 * t19;
t7 = -qJD(3) * pkin(3) + t40;
t6 = t19 * pkin(3) + t38;
t5 = -t19 * pkin(4) - t8;
t4 = t21 * pkin(4) - t50 * qJD(3) + t40;
t3 = t50 * t19 + t38;
t2 = t49 * t3 + t33 * t4;
t1 = -t33 * t3 + t49 * t4;
t12 = [0, 0, 0, 0, 0, t51, 0, 0, 0, 0, t28 * t51, t31 * t36 * t32, 0, t29 * t51, 0, 0, -t27 * t45, t27 * t46, (t28 + t29) * t36 * qJ(2), t27 ^ 2 / 0.2e1 + (t29 / 0.2e1 + t28 / 0.2e1) * qJ(2) ^ 2 * t36, t41, -t48, t43, t42, -t44, t30, t9 * qJD(3) + t25 * t19, -t10 * qJD(3) + t25 * t21, -t10 * t19 - t9 * t21, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t30, -t43, t44, t41, -t48, t42, t8 * t19 + t7 * t21, t7 * qJD(3) - t6 * t19, -t8 * qJD(3) - t6 * t21, t6 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t15, t11 ^ 2 / 0.2e1, -t11 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t5 * t11, t5 * t13 - t2 * t15, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t12;
