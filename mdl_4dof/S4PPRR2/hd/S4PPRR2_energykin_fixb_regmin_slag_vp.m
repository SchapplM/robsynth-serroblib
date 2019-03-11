% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% 
% Output:
% T_reg [1x8]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:00
% EndTime: 2019-03-08 18:17:00
% DurationCPUTime: 0.02s
% Computational Cost: add. (17->10), mult. (50->28), div. (0->0), fcn. (30->6), ass. (0->12)
t34 = sin(pkin(6));
t35 = cos(pkin(6));
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t41 = (-t34 * t37 + t35 * t39) * qJD(1);
t40 = qJD(1) ^ 2;
t38 = cos(qJ(4));
t36 = sin(qJ(4));
t33 = qJD(3) + qJD(4);
t31 = (t34 * t39 + t35 * t37) * qJD(1);
t30 = qJD(3) * pkin(3) + t41;
t1 = [t40 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t34 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1) * t40, qJD(3) ^ 2 / 0.2e1, t41 * qJD(3), -t31 * qJD(3), t33 ^ 2 / 0.2e1 (t38 * t30 - t36 * t31) * t33 -(t36 * t30 + t38 * t31) * t33;];
T_reg  = t1;
