% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:03
% EndTime: 2019-03-08 18:36:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (52->13), mult. (98->33), div. (0->0), fcn. (36->4), ass. (0->18)
t10 = sin(qJ(3));
t12 = cos(qJ(3));
t11 = sin(qJ(2));
t18 = pkin(1) * qJD(1);
t17 = t11 * t18;
t13 = cos(qJ(2));
t16 = t13 * t18;
t9 = qJD(1) + qJD(2);
t6 = t9 * pkin(2) + t16;
t4 = t10 * t6 + t12 * t17;
t8 = qJD(3) + t9;
t19 = t4 * t8;
t3 = -t10 * t17 + t12 * t6;
t14 = qJD(1) ^ 2;
t7 = t8 ^ 2 / 0.2e1;
t2 = t4 ^ 2 / 0.2e1;
t1 = t8 * pkin(3) + t3;
t5 = [0, 0, 0, 0, 0, t14 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 ^ 2 / 0.2e1, t9 * t16, -t9 * t17, 0 (t11 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t14, 0, 0, 0, 0, 0, t7, t3 * t8, -t19, 0, t2 + t3 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t7, t1 * t8, -t19, 0, t2 + t1 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t5;
