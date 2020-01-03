% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR11
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:53
% EndTime: 2019-12-31 18:27:54
% DurationCPUTime: 0.19s
% Computational Cost: add. (260->46), mult. (723->102), div. (0->0), fcn. (463->6), ass. (0->42)
t34 = sin(pkin(8));
t35 = cos(pkin(8));
t37 = sin(qJ(3));
t50 = cos(qJ(3));
t22 = (t50 * t34 + t35 * t37) * qJD(1);
t25 = -(pkin(2) * t35 + pkin(1)) * qJD(1) + qJD(2);
t53 = -t22 * qJ(4) + t25;
t38 = qJD(1) ^ 2;
t52 = t38 / 0.2e1;
t51 = -pkin(3) - pkin(4);
t49 = cos(qJ(5));
t44 = qJD(1) * t35;
t45 = qJD(1) * t34;
t20 = t37 * t45 - t50 * t44;
t48 = t22 * t20;
t47 = pkin(6) + qJ(2);
t23 = t47 * t45;
t24 = t47 * t44;
t13 = -t37 * t23 + t50 * t24;
t43 = qJD(3) * t20;
t42 = t20 ^ 2 / 0.2e1;
t11 = qJD(3) * qJ(4) + t13;
t12 = -t50 * t23 - t37 * t24;
t40 = qJD(4) - t12;
t36 = sin(qJ(5));
t32 = qJD(3) ^ 2 / 0.2e1;
t30 = qJD(3) - qJD(5);
t29 = t35 ^ 2;
t28 = t34 ^ 2;
t27 = -qJD(1) * pkin(1) + qJD(2);
t16 = t22 * qJD(3);
t15 = t22 ^ 2 / 0.2e1;
t10 = -qJD(3) * pkin(3) + t40;
t9 = t36 * t20 + t49 * t22;
t7 = -t49 * t20 + t36 * t22;
t6 = t20 * pkin(3) + t53;
t5 = t20 * pkin(7) + t11;
t4 = -t22 * pkin(7) + t51 * qJD(3) + t40;
t3 = t51 * t20 - t53;
t2 = t36 * t4 + t49 * t5;
t1 = -t36 * t5 + t49 * t4;
t8 = [0, 0, 0, 0, 0, t52, 0, 0, 0, 0, t28 * t52, t34 * t38 * t35, 0, t29 * t52, 0, 0, -t27 * t44, t27 * t45, (t28 + t29) * t38 * qJ(2), t27 ^ 2 / 0.2e1 + (t29 / 0.2e1 + t28 / 0.2e1) * qJ(2) ^ 2 * t38, t15, -t48, t16, t42, -t43, t32, t12 * qJD(3) + t25 * t20, -t13 * qJD(3) + t25 * t22, -t12 * t22 - t13 * t20, t13 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1, t15, t16, t48, t32, t43, t42, -t10 * qJD(3) + t6 * t20, t10 * t22 - t11 * t20, t11 * qJD(3) - t6 * t22, t11 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, -t9 * t30, t7 ^ 2 / 0.2e1, t7 * t30, t30 ^ 2 / 0.2e1, -t1 * t30 + t3 * t7, t2 * t30 + t3 * t9, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t8;
