% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:23
% EndTime: 2019-07-18 13:27:24
% DurationCPUTime: 0.02s
% Computational Cost: add. (17->7), mult. (30->20), div. (0->0), fcn. (10->4), ass. (0->9)
t36 = pkin(1) * qJD(2);
t29 = qJD(2) + qJD(3);
t35 = sin(qJ(3)) * t36;
t34 = cos(qJ(3)) * t36;
t32 = cos(qJ(4));
t30 = sin(qJ(4));
t28 = qJD(4) + t29;
t27 = t29 * pkin(2) + t34;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t29 ^ 2 / 0.2e1, t29 * t34, -t29 * t35, t28 ^ 2 / 0.2e1, (t32 * t27 - t30 * t35) * t28, -(t30 * t27 + t32 * t35) * t28;];
T_reg  = t1;
