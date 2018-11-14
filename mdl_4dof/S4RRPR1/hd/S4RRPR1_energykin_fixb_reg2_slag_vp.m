% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:29
% EndTime: 2018-11-14 13:53:29
% DurationCPUTime: 0.07s
% Computational Cost: add. (68->16), mult. (138->40), div. (0->0), fcn. (60->6), ass. (0->21)
t22 = pkin(1) * qJD(1);
t10 = qJD(1) + qJD(2);
t15 = sin(qJ(2));
t21 = t15 * t22;
t17 = cos(qJ(2));
t20 = t17 * t22;
t12 = sin(pkin(7));
t13 = cos(pkin(7));
t7 = pkin(2) * t10 + t20;
t4 = -t12 * t21 + t13 * t7;
t18 = qJD(1) ^ 2;
t16 = cos(qJ(4));
t14 = sin(qJ(4));
t11 = qJD(3) ^ 2 / 0.2e1;
t9 = qJD(4) + t10;
t8 = t10 ^ 2 / 0.2e1;
t5 = t12 * t7 + t13 * t21;
t3 = pkin(3) * t10 + t4;
t2 = t14 * t3 + t16 * t5;
t1 = -t14 * t5 + t16 * t3;
t6 = [0, 0, 0, 0, 0, t18 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10 * t20, -t10 * t21, 0 (t15 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t18, 0, 0, 0, 0, 0, t8, t4 * t10, -t5 * t10, 0, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t11, 0, 0, 0, 0, 0, t9 ^ 2 / 0.2e1, t1 * t9, -t2 * t9, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t11;];
T_reg  = t6;
