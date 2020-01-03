% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:45
% EndTime: 2019-12-31 20:55:45
% DurationCPUTime: 0.13s
% Computational Cost: add. (259->36), mult. (628->79), div. (0->0), fcn. (411->6), ass. (0->34)
t142 = -pkin(7) - pkin(6);
t131 = qJD(1) ^ 2;
t141 = t131 / 0.2e1;
t140 = cos(qJ(3));
t130 = cos(qJ(2));
t139 = t130 * t131;
t128 = sin(qJ(3));
t129 = sin(qJ(2));
t119 = (t128 * t130 + t140 * t129) * qJD(1);
t125 = qJD(2) + qJD(3);
t137 = qJD(1) * t129;
t121 = qJD(2) * pkin(2) + t142 * t137;
t136 = qJD(1) * t130;
t122 = t142 * t136;
t132 = t140 * t121 + t128 * t122;
t109 = t125 * pkin(3) - t119 * qJ(4) + t132;
t118 = t128 * t137 - t140 * t136;
t138 = t128 * t121 - t140 * t122;
t111 = -t118 * qJ(4) + t138;
t126 = sin(pkin(8));
t127 = cos(pkin(8));
t106 = t126 * t109 + t127 * t111;
t135 = qJD(1) * qJD(2);
t134 = t129 * t135;
t133 = t130 * t135;
t123 = (-pkin(2) * t130 - pkin(1)) * qJD(1);
t105 = t127 * t109 - t126 * t111;
t114 = t118 * pkin(3) + qJD(4) + t123;
t113 = -t126 * t118 + t127 * t119;
t112 = t127 * t118 + t126 * t119;
t107 = t112 * pkin(4) - t113 * qJ(5) + t114;
t104 = t125 * qJ(5) + t106;
t103 = -t125 * pkin(4) + qJD(5) - t105;
t1 = [t141, 0, 0, t129 ^ 2 * t141, t129 * t139, t134, t133, qJD(2) ^ 2 / 0.2e1, pkin(1) * t139 - pkin(6) * t134, -t131 * pkin(1) * t129 - pkin(6) * t133, t119 ^ 2 / 0.2e1, -t119 * t118, t119 * t125, -t118 * t125, t125 ^ 2 / 0.2e1, t123 * t118 + t132 * t125, t123 * t119 - t138 * t125, -t105 * t113 - t106 * t112, t106 ^ 2 / 0.2e1 + t105 ^ 2 / 0.2e1 + t114 ^ 2 / 0.2e1, -t103 * t125 + t107 * t112, t103 * t113 - t104 * t112, t104 * t125 - t107 * t113, t104 ^ 2 / 0.2e1 + t107 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1;];
T_reg = t1;
