% Calculate minimal parameter regressor of fixed base kinetic energy for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S2RR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_energykin_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_energykin_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_energykin_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:48:47
% EndTime: 2018-11-16 16:48:47
% DurationCPUTime: 0.02s
% Computational Cost: add. (2->2), mult. (21->13), div. (0->0), fcn. (7->2), ass. (0->8)
t20 = qJD(1) ^ 2;
t24 = t20 / 0.2e1;
t23 = qJD(1) * qJD(2);
t18 = sin(qJ(2));
t22 = t18 * t23;
t19 = cos(qJ(2));
t21 = t19 * t23;
t1 = [t24, 0, 0, t18 ^ 2 * t24, t18 * t20 * t19, t22, t21, qJD(2) ^ 2 / 0.2e1, -pkin(1) * t22, -pkin(1) * t21;];
T_reg  = t1;
