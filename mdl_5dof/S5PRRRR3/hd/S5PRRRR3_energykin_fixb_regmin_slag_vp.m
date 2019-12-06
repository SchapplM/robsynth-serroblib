% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:18
% EndTime: 2019-12-05 17:06:18
% DurationCPUTime: 0.06s
% Computational Cost: add. (57->13), mult. (85->40), div. (0->0), fcn. (37->6), ass. (0->18)
t67 = qJD(2) + qJD(3);
t66 = qJD(4) + t67;
t65 = t66 ^ 2;
t81 = t65 / 0.2e1;
t78 = pkin(2) * qJD(2);
t75 = cos(qJ(3)) * t78;
t63 = t67 * pkin(3) + t75;
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t76 = sin(qJ(3)) * t78;
t74 = t72 * t63 - t69 * t76;
t80 = (-t66 * pkin(4) - t74) * t66;
t79 = t69 * t63 + t72 * t76;
t77 = qJD(5) * t66;
t71 = cos(qJ(5));
t68 = sin(qJ(5));
t61 = t66 * pkin(8) + t79;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t67 ^ 2 / 0.2e1, t67 * t75, -t67 * t76, t81, t74 * t66, -t79 * t66, t68 ^ 2 * t81, t68 * t65 * t71, t68 * t77, t71 * t77, qJD(5) ^ 2 / 0.2e1, -t71 * t80 + (t71 * qJD(1) - t68 * t61) * qJD(5), t68 * t80 - (t68 * qJD(1) + t71 * t61) * qJD(5);];
T_reg = t1;
