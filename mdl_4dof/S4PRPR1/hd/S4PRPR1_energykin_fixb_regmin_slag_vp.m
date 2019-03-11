% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:02
% EndTime: 2019-03-08 18:21:02
% DurationCPUTime: 0.02s
% Computational Cost: add. (15->10), mult. (31->20), div. (0->0), fcn. (4->2), ass. (0->10)
t39 = qJD(2) ^ 2;
t41 = t39 / 0.2e1;
t40 = qJ(3) * qJD(2);
t38 = cos(qJ(4));
t37 = sin(qJ(4));
t36 = qJD(1) ^ 2 / 0.2e1;
t35 = -qJD(2) + qJD(4);
t34 = -qJD(2) * pkin(2) + qJD(3);
t33 = qJD(3) + (-pkin(2) - pkin(3)) * qJD(2);
t1 = [t36, t41, 0, 0, -t34 * qJD(2), t39 * qJ(3), qJ(3) ^ 2 * t41 + t36 + t34 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1 (t38 * t33 - t37 * t40) * t35 -(t37 * t33 + t38 * t40) * t35;];
T_reg  = t1;
