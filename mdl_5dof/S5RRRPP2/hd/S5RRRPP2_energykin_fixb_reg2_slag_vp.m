% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:57
% EndTime: 2019-12-31 20:51:58
% DurationCPUTime: 0.14s
% Computational Cost: add. (171->37), mult. (279->81), div. (0->0), fcn. (98->4), ass. (0->37)
t22 = sin(qJ(3));
t20 = t22 ^ 2;
t41 = t20 / 0.2e1;
t24 = cos(qJ(3));
t21 = t24 ^ 2;
t40 = t21 / 0.2e1;
t39 = pkin(3) + pkin(4);
t17 = qJD(1) + qJD(2);
t23 = sin(qJ(2));
t36 = pkin(1) * qJD(1);
t30 = t23 * t36;
t10 = t17 * pkin(7) + t30;
t6 = qJD(3) * qJ(4) + t24 * t10;
t38 = t17 * t22;
t37 = t17 * t24;
t35 = qJ(5) * t17;
t34 = t22 * t10 + qJD(4);
t33 = qJD(3) * t22;
t32 = qJD(3) * t24;
t16 = t17 ^ 2;
t31 = t22 * t16 * t24;
t25 = cos(qJ(2));
t29 = t25 * t36;
t28 = qJ(4) * t22 + pkin(2);
t26 = qJD(1) ^ 2;
t18 = qJD(3) ^ 2 / 0.2e1;
t15 = t17 * t32;
t14 = t17 * t33;
t13 = t16 * t40;
t12 = t16 * t41;
t11 = -t17 * pkin(2) - t29;
t5 = -qJD(3) * pkin(3) + t34;
t4 = -t24 * t35 + t6;
t3 = -t29 + (-pkin(3) * t24 - t28) * t17;
t2 = -t39 * qJD(3) - t22 * t35 + t34;
t1 = t29 + qJD(5) + (t39 * t24 + t28) * t17;
t7 = [0, 0, 0, 0, 0, t26 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 / 0.2e1, t17 * t29, -t17 * t30, 0, (t23 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t26, t12, t31, t14, t13, t15, t18, -t10 * t33 - t11 * t37, -t10 * t32 + t11 * t38, (t20 + t21) * t17 * t10, t11 ^ 2 / 0.2e1 + (t40 + t41) * t10 ^ 2, t12, t14, -t31, t18, -t15, t13, -t5 * qJD(3) - t3 * t37, (t22 * t5 + t24 * t6) * t17, t6 * qJD(3) - t3 * t38, t6 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t12, -t31, -t14, t13, t15, t18, -t2 * qJD(3) + t1 * t37, t4 * qJD(3) + t1 * t38, (-t2 * t22 - t24 * t4) * t17, t4 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg = t7;
