% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PPPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:10:17
% EndTime: 2019-03-08 18:10:17
% DurationCPUTime: 0.05s
% Computational Cost: add. (16->9), mult. (47->22), div. (0->0), fcn. (20->4), ass. (0->12)
t13 = qJD(1) ^ 2 / 0.2e1;
t5 = qJD(2) ^ 2 / 0.2e1;
t6 = sin(pkin(5));
t12 = t6 ^ 2 * t13 + t5;
t11 = qJD(1) * t6;
t9 = cos(qJ(4));
t8 = sin(qJ(4));
t7 = cos(pkin(5));
t3 = -t7 * qJD(1) + qJD(3);
t2 = t9 * t11 + t8 * t3;
t1 = -t8 * t11 + t9 * t3;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 * t13 + t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 / 0.2e1 + t12, 0, 0, 0, 0, 0, qJD(4) ^ 2 / 0.2e1, t1 * qJD(4), -t2 * qJD(4), 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5;];
T_reg  = t4;
