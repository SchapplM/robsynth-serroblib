% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% 
% Output:
% T_reg [1x8]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:24
% EndTime: 2019-03-08 18:12:24
% DurationCPUTime: 0.02s
% Computational Cost: add. (11->8), mult. (28->18), div. (0->0), fcn. (6->2), ass. (0->7)
t39 = qJD(2) * qJD(3);
t38 = cos(qJ(3));
t37 = sin(qJ(3));
t36 = qJD(1) ^ 2 / 0.2e1;
t35 = qJD(3) * qJ(4) + t37 * qJD(2);
t34 = -qJD(3) * pkin(3) - t38 * qJD(2) + qJD(4);
t1 = [t36, t36 + qJD(2) ^ 2 / 0.2e1, qJD(3) ^ 2 / 0.2e1, t38 * t39, -t37 * t39, -t34 * qJD(3), t35 * qJD(3), t35 ^ 2 / 0.2e1 + t36 + t34 ^ 2 / 0.2e1;];
T_reg  = t1;
