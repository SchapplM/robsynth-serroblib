% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:31
% EndTime: 2019-12-05 16:13:31
% DurationCPUTime: 0.07s
% Computational Cost: add. (136->30), mult. (293->70), div. (0->0), fcn. (163->6), ass. (0->27)
t132 = qJD(2) ^ 2;
t139 = t132 / 0.2e1;
t128 = sin(qJ(3));
t130 = cos(qJ(3));
t131 = cos(qJ(2));
t135 = t131 * qJD(1);
t116 = -t135 + (-pkin(3) * t130 - qJ(4) * t128 - pkin(2)) * qJD(2);
t129 = sin(qJ(2));
t123 = qJD(2) * pkin(6) + t129 * qJD(1);
t119 = qJD(3) * qJ(4) + t130 * t123;
t126 = sin(pkin(8));
t127 = cos(pkin(8));
t114 = t126 * t116 + t127 * t119;
t138 = qJD(2) * t128;
t137 = qJD(2) * t130;
t136 = qJD(3) * t123;
t134 = qJD(1) * qJD(2);
t133 = qJD(2) * qJD(3);
t118 = -qJD(3) * pkin(3) + t128 * t123 + qJD(4);
t113 = t127 * t116 - t126 * t119;
t124 = -qJD(2) * pkin(2) - t135;
t121 = t126 * qJD(3) + t127 * t138;
t120 = -t127 * qJD(3) + t126 * t138;
t112 = t120 * pkin(4) - t121 * qJ(5) + t118;
t111 = -qJ(5) * t137 + t114;
t110 = pkin(4) * t137 + qJD(5) - t113;
t1 = [qJD(1) ^ 2 / 0.2e1, t139, t131 * t134, -t129 * t134, t128 ^ 2 * t139, t128 * t132 * t130, t128 * t133, t130 * t133, qJD(3) ^ 2 / 0.2e1, -t124 * t137 - t128 * t136, t124 * t138 - t130 * t136, -t113 * t137 + t118 * t120, t114 * t137 + t118 * t121, -t113 * t121 - t114 * t120, t114 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1, t110 * t137 + t112 * t120, t110 * t121 - t111 * t120, -t111 * t137 - t112 * t121, t111 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1 + t110 ^ 2 / 0.2e1;];
T_reg = t1;
