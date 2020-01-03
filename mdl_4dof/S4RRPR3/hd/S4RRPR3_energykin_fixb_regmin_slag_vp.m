% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR3
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
% T_reg [1x14]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:33
% EndTime: 2019-12-31 17:01:33
% DurationCPUTime: 0.04s
% Computational Cost: add. (42->13), mult. (85->40), div. (0->0), fcn. (37->6), ass. (0->17)
t61 = qJD(1) + qJD(2);
t60 = t61 ^ 2;
t73 = t60 / 0.2e1;
t71 = pkin(1) * qJD(1);
t68 = cos(qJ(2)) * t71;
t58 = t61 * pkin(2) + t68;
t62 = sin(pkin(7));
t63 = cos(pkin(7));
t69 = sin(qJ(2)) * t71;
t55 = t63 * t58 - t62 * t69;
t72 = (-t61 * pkin(3) - t55) * t61;
t56 = t62 * t58 + t63 * t69;
t70 = qJD(4) * t61;
t66 = cos(qJ(4));
t64 = sin(qJ(4));
t54 = t61 * pkin(6) + t56;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t73, t61 * t68, -t61 * t69, t56 ^ 2 / 0.2e1 + t55 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t64 ^ 2 * t73, t64 * t60 * t66, t64 * t70, t66 * t70, qJD(4) ^ 2 / 0.2e1, -t66 * t72 + (t66 * qJD(3) - t64 * t54) * qJD(4), t64 * t72 - (t64 * qJD(3) + t66 * t54) * qJD(4);];
T_reg = t1;
