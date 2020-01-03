% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:11
% EndTime: 2019-12-31 18:11:11
% DurationCPUTime: 0.11s
% Computational Cost: add. (112->37), mult. (297->81), div. (0->0), fcn. (118->4), ass. (0->32)
t25 = qJD(1) ^ 2;
t19 = t25 / 0.2e1;
t36 = pkin(3) + pkin(4);
t35 = pkin(1) * t25;
t21 = sin(pkin(7));
t11 = (pkin(1) * t21 + pkin(6)) * qJD(1);
t23 = sin(qJ(3));
t24 = cos(qJ(3));
t8 = t23 * qJD(2) + t24 * t11;
t34 = qJD(1) * t23;
t33 = qJD(1) * t24;
t32 = qJ(5) * qJD(1);
t31 = qJD(1) * qJD(3);
t30 = t23 * t25 * t24;
t5 = qJD(3) * qJ(4) + t8;
t22 = cos(pkin(7));
t29 = -pkin(1) * t22 - pkin(2);
t7 = t24 * qJD(2) - t23 * t11;
t28 = qJD(4) - t7;
t27 = qJ(4) * t23 - t29;
t18 = qJD(3) ^ 2 / 0.2e1;
t16 = t24 * t31;
t15 = t23 * t31;
t14 = t24 ^ 2 * t19;
t13 = t23 ^ 2 * t19;
t12 = t29 * qJD(1);
t6 = (-pkin(3) * t24 - t27) * qJD(1);
t4 = -qJD(3) * pkin(3) + t28;
t3 = -t24 * t32 + t5;
t2 = qJD(5) + (t36 * t24 + t27) * qJD(1);
t1 = -t36 * qJD(3) - t23 * t32 + t28;
t9 = [0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t22 * t35, -t21 * t35, 0, qJD(2) ^ 2 / 0.2e1 + (t21 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t25, t13, t30, t15, t14, t16, t18, t7 * qJD(3) - t12 * t33, -t8 * qJD(3) + t12 * t34, (-t23 * t7 + t24 * t8) * qJD(1), t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t13, t15, -t30, t18, -t16, t14, -t4 * qJD(3) - t6 * t33, (t23 * t4 + t24 * t5) * qJD(1), t5 * qJD(3) - t6 * t34, t5 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t13, -t30, -t15, t14, t16, t18, -t1 * qJD(3) + t2 * t33, t3 * qJD(3) + t2 * t34, (-t1 * t23 - t24 * t3) * qJD(1), t3 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg = t9;
