% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:56
% EndTime: 2019-12-31 19:54:56
% DurationCPUTime: 0.13s
% Computational Cost: add. (243->36), mult. (615->79), div. (0->0), fcn. (413->6), ass. (0->32)
t140 = qJD(1) * (pkin(6) + qJ(3));
t131 = qJD(1) ^ 2;
t139 = t131 / 0.2e1;
t130 = cos(qJ(2));
t137 = t130 * t131;
t128 = sin(qJ(2));
t121 = qJD(2) * pkin(2) - t128 * t140;
t122 = t130 * t140;
t125 = sin(pkin(8));
t126 = cos(pkin(8));
t112 = t126 * t121 - t125 * t122;
t119 = (t125 * t130 + t126 * t128) * qJD(1);
t108 = qJD(2) * pkin(3) - t119 * pkin(7) + t112;
t113 = t125 * t121 + t126 * t122;
t118 = (-t125 * t128 + t126 * t130) * qJD(1);
t109 = t118 * pkin(7) + t113;
t127 = sin(qJ(4));
t129 = cos(qJ(4));
t136 = t127 * t108 + t129 * t109;
t135 = qJD(1) * qJD(2);
t134 = t128 * t135;
t133 = t130 * t135;
t132 = t129 * t108 - t127 * t109;
t123 = qJD(3) + (-pkin(2) * t130 - pkin(1)) * qJD(1);
t114 = -t118 * pkin(3) + t123;
t124 = qJD(2) + qJD(4);
t111 = t127 * t118 + t129 * t119;
t110 = -t129 * t118 + t127 * t119;
t105 = t110 * pkin(4) - t111 * qJ(5) + t114;
t104 = t124 * qJ(5) + t136;
t103 = -t124 * pkin(4) + qJD(5) - t132;
t1 = [t139, 0, 0, t128 ^ 2 * t139, t128 * t137, t134, t133, qJD(2) ^ 2 / 0.2e1, pkin(1) * t137 - pkin(6) * t134, -t131 * pkin(1) * t128 - pkin(6) * t133, -t112 * t119 + t113 * t118, t113 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1, t111 ^ 2 / 0.2e1, -t111 * t110, t111 * t124, -t110 * t124, t124 ^ 2 / 0.2e1, t114 * t110 + t132 * t124, t114 * t111 - t136 * t124, -t103 * t124 + t105 * t110, t103 * t111 - t104 * t110, t104 * t124 - t105 * t111, t104 ^ 2 / 0.2e1 + t105 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1;];
T_reg = t1;
