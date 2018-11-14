% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:34
% EndTime: 2018-11-14 13:50:34
% DurationCPUTime: 0.07s
% Computational Cost: add. (54->17), mult. (140->41), div. (0->0), fcn. (60->6), ass. (0->20)
t18 = qJD(1) ^ 2;
t21 = pkin(1) * t18;
t9 = qJD(1) + qJD(3);
t12 = sin(pkin(7));
t20 = pkin(1) * qJD(1) * t12;
t15 = sin(qJ(3));
t17 = cos(qJ(3));
t13 = cos(pkin(7));
t7 = (pkin(1) * t13 + pkin(2)) * qJD(1);
t4 = -t15 * t20 + t17 * t7;
t16 = cos(qJ(4));
t14 = sin(qJ(4));
t11 = t18 / 0.2e1;
t10 = qJD(2) ^ 2 / 0.2e1;
t8 = qJD(4) + t9;
t5 = t15 * t7 + t17 * t20;
t3 = pkin(3) * t9 + t4;
t2 = t14 * t3 + t16 * t5;
t1 = -t14 * t5 + t16 * t3;
t6 = [0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t13 * t21, -t12 * t21, 0, t10 + (t12 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t18, 0, 0, 0, 0, 0, t9 ^ 2 / 0.2e1, t4 * t9, -t5 * t9, 0, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t10, 0, 0, 0, 0, 0, t8 ^ 2 / 0.2e1, t1 * t8, -t2 * t8, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10;];
T_reg  = t6;
