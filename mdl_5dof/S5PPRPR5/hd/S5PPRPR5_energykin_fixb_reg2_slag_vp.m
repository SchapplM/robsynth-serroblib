% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRPR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:29
% EndTime: 2019-12-31 17:33:29
% DurationCPUTime: 0.08s
% Computational Cost: add. (47->20), mult. (109->49), div. (0->0), fcn. (39->4), ass. (0->19)
t13 = qJD(3) ^ 2;
t7 = t13 / 0.2e1;
t10 = sin(qJ(3));
t5 = qJD(3) * qJ(4) + qJD(2) * t10;
t19 = t5 * qJD(3);
t18 = t5 ^ 2 / 0.2e1;
t17 = qJD(2) * qJD(3);
t16 = qJD(3) * qJD(5);
t12 = cos(qJ(3));
t15 = -qJD(2) * t12 + qJD(4);
t14 = qJD(2) ^ 2;
t11 = cos(qJ(5));
t9 = sin(qJ(5));
t8 = qJD(1) ^ 2 / 0.2e1;
t4 = -qJD(3) * pkin(3) + t15;
t3 = (-pkin(3) - pkin(6)) * qJD(3) + t15;
t2 = -qJD(1) * t11 + t3 * t9;
t1 = qJD(1) * t9 + t11 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 + t14 / 0.2e1, 0, 0, 0, 0, 0, t7, t12 * t17, -t10 * t17, 0, t8 + (t10 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1) * t14, t7, 0, 0, 0, 0, 0, 0, t4 * qJD(3), t19, t8 + t18 + t4 ^ 2 / 0.2e1, t11 ^ 2 * t7, -t11 * t13 * t9, t11 * t16, t9 ^ 2 * t7, -t9 * t16, qJD(5) ^ 2 / 0.2e1, qJD(5) * t1 + t9 * t19, -qJD(5) * t2 + t11 * t19, (-t1 * t11 - t2 * t9) * qJD(3), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t18;];
T_reg = t6;
