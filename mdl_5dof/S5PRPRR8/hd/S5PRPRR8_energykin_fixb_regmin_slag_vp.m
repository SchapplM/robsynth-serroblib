% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:26
% EndTime: 2019-12-05 16:04:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (66->28), mult. (158->64), div. (0->0), fcn. (92->8), ass. (0->31)
t127 = qJD(2) ^ 2;
t141 = t127 / 0.2e1;
t140 = qJD(1) ^ 2 / 0.2e1;
t126 = cos(qJ(2));
t138 = qJD(1) * sin(pkin(5));
t130 = -t126 * t138 + qJD(3);
t110 = (-pkin(2) - pkin(7)) * qJD(2) + t130;
t122 = sin(qJ(4));
t125 = cos(qJ(4));
t120 = cos(pkin(5));
t137 = qJD(1) * t120;
t139 = t122 * t110 + t125 * t137;
t136 = qJD(2) * t125;
t123 = sin(qJ(2));
t132 = t123 * t138;
t114 = qJD(2) * qJ(3) + t132;
t135 = t114 * qJD(2);
t134 = t122 * qJD(2);
t133 = qJD(2) * qJD(4);
t131 = qJD(2) * t138;
t129 = t125 * t110 - t122 * t137;
t124 = cos(qJ(5));
t121 = sin(qJ(5));
t117 = qJD(5) + t134;
t113 = t121 * qJD(4) + t124 * t136;
t112 = -t124 * qJD(4) + t121 * t136;
t111 = -qJD(2) * pkin(2) + t130;
t108 = t132 + (pkin(4) * t122 - pkin(8) * t125 + qJ(3)) * qJD(2);
t107 = qJD(4) * pkin(8) + t139;
t106 = -qJD(4) * pkin(4) - t129;
t1 = [t140, t141, t126 * t131, -t123 * t131, t111 * qJD(2), t135, t120 ^ 2 * t140 + t114 ^ 2 / 0.2e1 + t111 ^ 2 / 0.2e1, t125 ^ 2 * t141, -t125 * t127 * t122, t125 * t133, -t122 * t133, qJD(4) ^ 2 / 0.2e1, t129 * qJD(4) + t114 * t134, -t139 * qJD(4) + t125 * t135, t113 ^ 2 / 0.2e1, -t113 * t112, t113 * t117, -t112 * t117, t117 ^ 2 / 0.2e1, (-t121 * t107 + t124 * t108) * t117 + t106 * t112, -(t124 * t107 + t121 * t108) * t117 + t106 * t113;];
T_reg = t1;
