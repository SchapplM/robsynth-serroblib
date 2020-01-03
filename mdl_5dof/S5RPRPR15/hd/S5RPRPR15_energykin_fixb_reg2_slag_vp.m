% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR15_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:23
% EndTime: 2019-12-31 18:37:23
% DurationCPUTime: 0.15s
% Computational Cost: add. (251->43), mult. (541->100), div. (0->0), fcn. (293->6), ass. (0->35)
t35 = qJD(1) ^ 2;
t28 = t35 / 0.2e1;
t42 = cos(qJ(5));
t33 = sin(qJ(3));
t34 = cos(qJ(3));
t16 = (pkin(3) * t33 - qJ(4) * t34 + qJ(2)) * qJD(1);
t22 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t17 = qJD(3) * qJ(4) + t33 * t22;
t31 = sin(pkin(8));
t40 = cos(pkin(8));
t6 = t31 * t16 + t40 * t17;
t41 = t35 * qJ(2);
t39 = qJD(1) * t34;
t38 = qJD(3) * t22;
t37 = t33 * qJD(1);
t36 = qJD(1) * qJD(3);
t5 = t40 * t16 - t31 * t17;
t15 = -qJD(3) * pkin(3) - t34 * t22 + qJD(4);
t32 = sin(qJ(5));
t30 = t34 ^ 2;
t29 = t33 ^ 2;
t27 = qJ(2) ^ 2 * t28;
t25 = -pkin(1) * qJD(1) + qJD(2);
t24 = t29 * t28;
t23 = qJD(5) + t37;
t20 = t31 * qJD(3) + t39 * t40;
t18 = -qJD(3) * t40 + t31 * t39;
t10 = t18 * pkin(4) + t15;
t9 = -t32 * t18 + t42 * t20;
t7 = t42 * t18 + t32 * t20;
t4 = -t18 * pkin(7) + t6;
t3 = pkin(4) * t37 - t20 * pkin(7) + t5;
t2 = t32 * t3 + t42 * t4;
t1 = t42 * t3 - t32 * t4;
t8 = [0, 0, 0, 0, 0, t28, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, t25 * qJD(1), t41, t27 + t25 ^ 2 / 0.2e1, t30 * t28, -t34 * t35 * t33, t34 * t36, t24, -t33 * t36, qJD(3) ^ 2 / 0.2e1, t33 * t41 + t34 * t38, -t33 * t38 + t34 * t41, (-t29 - t30) * t22 * qJD(1), t27 + (t29 / 0.2e1 + t30 / 0.2e1) * t22 ^ 2, t20 ^ 2 / 0.2e1, -t20 * t18, t20 * t37, t18 ^ 2 / 0.2e1, -t18 * t37, t24, t15 * t18 + t37 * t5, t15 * t20 - t37 * t6, -t6 * t18 - t5 * t20, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t23, t7 ^ 2 / 0.2e1, -t7 * t23, t23 ^ 2 / 0.2e1, t1 * t23 + t10 * t7, t10 * t9 - t2 * t23, -t1 * t9 - t2 * t7, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t8;
