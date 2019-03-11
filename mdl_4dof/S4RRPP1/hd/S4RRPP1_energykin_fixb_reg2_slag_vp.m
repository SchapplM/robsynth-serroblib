% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:09
% EndTime: 2019-03-08 18:33:09
% DurationCPUTime: 0.06s
% Computational Cost: add. (51->15), mult. (102->35), div. (0->0), fcn. (36->4), ass. (0->17)
t11 = sin(pkin(6));
t12 = cos(pkin(6));
t13 = sin(qJ(2));
t19 = pkin(1) * qJD(1);
t18 = t13 * t19;
t14 = cos(qJ(2));
t17 = t14 * t19;
t9 = qJD(1) + qJD(2);
t6 = t9 * pkin(2) + t17;
t4 = t11 * t6 + t12 * t18;
t3 = -t11 * t18 + t12 * t6;
t15 = qJD(1) ^ 2;
t10 = qJD(3) ^ 2 / 0.2e1;
t8 = t9 ^ 2 / 0.2e1;
t2 = t9 * qJ(4) + t4;
t1 = -t9 * pkin(3) + qJD(4) - t3;
t5 = [0, 0, 0, 0, 0, t15 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9 * t17, -t9 * t18, 0 (t13 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t15, 0, 0, 0, 0, 0, t8, t3 * t9, -t4 * t9, 0, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t10, 0, 0, 0, t8, 0, 0, -t1 * t9, 0, t2 * t9, t2 ^ 2 / 0.2e1 + t10 + t1 ^ 2 / 0.2e1;];
T_reg  = t5;
