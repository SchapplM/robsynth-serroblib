% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:50
% EndTime: 2019-12-05 17:04:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (51->13), mult. (83->38), div. (0->0), fcn. (37->6), ass. (0->17)
t63 = qJD(2) + qJD(3);
t62 = qJD(4) + t63;
t61 = t62 ^ 2;
t76 = t61 / 0.2e1;
t73 = pkin(2) * qJD(2);
t70 = cos(qJ(3)) * t73;
t59 = t63 * pkin(3) + t70;
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t71 = sin(qJ(3)) * t73;
t75 = (-t68 * t59 + t65 * t71) * t62;
t74 = t65 * t59 + t68 * t71;
t72 = qJD(5) * t62;
t67 = cos(qJ(5));
t64 = sin(qJ(5));
t56 = t62 * pkin(6) + t74;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t63 ^ 2 / 0.2e1, t63 * t70, -t63 * t71, t76, -t75, -t74 * t62, t64 ^ 2 * t76, t64 * t61 * t67, t64 * t72, t67 * t72, qJD(5) ^ 2 / 0.2e1, -t67 * t75 + (t67 * qJD(1) - t64 * t56) * qJD(5), t64 * t75 - (t64 * qJD(1) + t67 * t56) * qJD(5);];
T_reg = t1;
