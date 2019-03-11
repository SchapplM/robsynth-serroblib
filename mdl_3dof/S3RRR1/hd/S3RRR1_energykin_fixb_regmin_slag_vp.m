% Calculate minimal parameter regressor of fixed base kinetic energy for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% T_reg [1x9]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S3RRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energykin_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_energykin_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:07
% EndTime: 2019-03-08 18:08:07
% DurationCPUTime: 0.02s
% Computational Cost: add. (17->7), mult. (28->18), div. (0->0), fcn. (10->4), ass. (0->9)
t35 = pkin(1) * qJD(1);
t28 = qJD(1) + qJD(2);
t34 = sin(qJ(2)) * t35;
t33 = cos(qJ(2)) * t35;
t31 = cos(qJ(3));
t29 = sin(qJ(3));
t27 = qJD(3) + t28;
t26 = t28 * pkin(2) + t33;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t28 ^ 2 / 0.2e1, t28 * t33, -t28 * t34, t27 ^ 2 / 0.2e1 (t31 * t26 - t29 * t34) * t27 -(t29 * t26 + t31 * t34) * t27;];
T_reg  = t1;
