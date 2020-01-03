% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x14]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:38
% EndTime: 2019-12-31 16:31:38
% DurationCPUTime: 0.03s
% Computational Cost: add. (22->8), mult. (47->31), div. (0->0), fcn. (17->4), ass. (0->12)
t46 = qJD(2) + qJD(3);
t45 = t46 ^ 2;
t56 = t45 / 0.2e1;
t54 = pkin(2) * qJD(2);
t51 = cos(qJ(3)) * t54;
t55 = (-t46 * pkin(3) - t51) * t46;
t53 = qJD(4) * t46;
t52 = sin(qJ(3)) * t54;
t49 = cos(qJ(4));
t47 = sin(qJ(4));
t43 = t46 * pkin(6) + t52;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t56, t46 * t51, -t46 * t52, t47 ^ 2 * t56, t47 * t45 * t49, t47 * t53, t49 * t53, qJD(4) ^ 2 / 0.2e1, -t49 * t55 + (t49 * qJD(1) - t47 * t43) * qJD(4), t47 * t55 - (t47 * qJD(1) + t49 * t43) * qJD(4);];
T_reg = t1;
