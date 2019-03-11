% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:59
% EndTime: 2019-03-08 18:22:59
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->10), mult. (50->24), div. (0->0), fcn. (8->2), ass. (0->12)
t12 = pkin(2) * qJD(2);
t6 = sin(qJ(3));
t11 = t6 * t12;
t7 = cos(qJ(3));
t10 = t7 * t12;
t8 = qJD(2) ^ 2;
t5 = qJD(1) ^ 2 / 0.2e1;
t4 = qJD(2) + qJD(3);
t3 = t4 ^ 2 / 0.2e1;
t2 = qJ(4) * t4 + t11;
t1 = -pkin(3) * t4 + qJD(4) - t10;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, t8 / 0.2e1, 0, 0, 0, t5, 0, 0, 0, 0, 0, t3, t4 * t10, -t4 * t11, 0, t5 + (t6 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t8, 0, 0, 0, t3, 0, 0, -t1 * t4, 0, t2 * t4, t2 ^ 2 / 0.2e1 + t5 + t1 ^ 2 / 0.2e1;];
T_reg  = t9;
