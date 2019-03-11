% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:00
% EndTime: 2019-03-08 18:35:00
% DurationCPUTime: 0.02s
% Computational Cost: add. (37->12), mult. (68->29), div. (0->0), fcn. (30->6), ass. (0->14)
t62 = pkin(1) * qJD(1);
t53 = qJD(1) + qJD(2);
t61 = sin(qJ(2)) * t62;
t60 = cos(qJ(2)) * t62;
t51 = t53 * pkin(2) + t60;
t54 = sin(pkin(7));
t55 = cos(pkin(7));
t48 = t55 * t51 - t54 * t61;
t58 = cos(qJ(4));
t56 = sin(qJ(4));
t52 = qJD(4) + t53;
t49 = t54 * t51 + t55 * t61;
t47 = t53 * pkin(3) + t48;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t53 ^ 2 / 0.2e1, t53 * t60, -t53 * t61, t49 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t52 ^ 2 / 0.2e1 (t58 * t47 - t56 * t49) * t52 -(t56 * t47 + t58 * t49) * t52;];
T_reg  = t1;
