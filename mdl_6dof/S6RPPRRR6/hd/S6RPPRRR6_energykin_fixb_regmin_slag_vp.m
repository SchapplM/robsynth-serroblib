% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:28
% EndTime: 2019-03-09 02:31:28
% DurationCPUTime: 0.07s
% Computational Cost: add. (148->38), mult. (283->83), div. (0->0), fcn. (154->6), ass. (0->33)
t135 = qJD(1) ^ 2;
t145 = t135 / 0.2e1;
t144 = cos(qJ(5));
t143 = -pkin(1) - qJ(3);
t132 = sin(qJ(4));
t134 = cos(qJ(4));
t114 = -qJD(2) + (pkin(4) * t132 - pkin(8) * t134 - t143) * qJD(1);
t126 = qJD(1) * qJ(2) + qJD(3);
t121 = -qJD(1) * pkin(7) + t126;
t117 = qJD(4) * pkin(8) + t132 * t121;
t131 = sin(qJ(5));
t142 = t131 * t114 + t144 * t117;
t141 = qJD(1) * t134;
t140 = qJD(4) * t121;
t122 = -t143 * qJD(1) - qJD(2);
t139 = t122 * qJD(1);
t138 = t132 * qJD(1);
t137 = qJD(1) * qJD(4);
t136 = t144 * t114 - t131 * t117;
t125 = qJD(5) + t138;
t118 = -qJD(4) * pkin(4) - t134 * t121;
t133 = cos(qJ(6));
t130 = sin(qJ(6));
t127 = -qJD(1) * pkin(1) + qJD(2);
t124 = qJD(6) + t125;
t120 = t131 * qJD(4) + t144 * t141;
t119 = -t144 * qJD(4) + t131 * t141;
t111 = t119 * pkin(5) + t118;
t110 = -t130 * t119 + t133 * t120;
t109 = t133 * t119 + t130 * t120;
t108 = -t119 * pkin(9) + t142;
t107 = t125 * pkin(5) - t120 * pkin(9) + t136;
t1 = [t145, 0, 0, t127 * qJD(1), t135 * qJ(2), qJ(2) ^ 2 * t145 + t127 ^ 2 / 0.2e1, t126 * qJD(1), t139, t122 ^ 2 / 0.2e1 + t126 ^ 2 / 0.2e1, t134 ^ 2 * t145, -t134 * t135 * t132, t134 * t137, -t132 * t137, qJD(4) ^ 2 / 0.2e1, t122 * t138 + t134 * t140, -t132 * t140 + t134 * t139, t120 ^ 2 / 0.2e1, -t120 * t119, t120 * t125, -t119 * t125, t125 ^ 2 / 0.2e1, t118 * t119 + t136 * t125, t118 * t120 - t142 * t125, t110 ^ 2 / 0.2e1, -t110 * t109, t110 * t124, -t109 * t124, t124 ^ 2 / 0.2e1 (t133 * t107 - t130 * t108) * t124 + t111 * t109 -(t130 * t107 + t133 * t108) * t124 + t111 * t110;];
T_reg  = t1;
