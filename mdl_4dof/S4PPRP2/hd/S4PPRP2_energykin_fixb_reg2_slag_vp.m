% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PPRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:13:20
% EndTime: 2019-03-08 18:13:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->13), mult. (68->29), div. (0->0), fcn. (34->4), ass. (0->12)
t10 = cos(pkin(5));
t11 = sin(qJ(3));
t12 = cos(qJ(3));
t9 = sin(pkin(5));
t3 = (t10 * t12 - t11 * t9) * qJD(1);
t4 = (t10 * t11 + t12 * t9) * qJD(1);
t13 = qJD(1) ^ 2;
t8 = qJD(2) ^ 2 / 0.2e1;
t7 = qJD(3) ^ 2 / 0.2e1;
t2 = qJD(3) * qJ(4) + t4;
t1 = -qJD(3) * pkin(3) + qJD(4) - t3;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t13 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 + (t9 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1) * t13, 0, 0, 0, 0, 0, t7, t3 * qJD(3), -t4 * qJD(3), 0, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t8, 0, 0, 0, t7, 0, 0, -t1 * qJD(3), 0, t2 * qJD(3), t2 ^ 2 / 0.2e1 + t8 + t1 ^ 2 / 0.2e1;];
T_reg  = t5;
