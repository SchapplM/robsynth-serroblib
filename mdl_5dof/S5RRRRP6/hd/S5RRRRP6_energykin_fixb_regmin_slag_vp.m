% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:33
% EndTime: 2019-12-31 21:54:33
% DurationCPUTime: 0.08s
% Computational Cost: add. (185->35), mult. (430->78), div. (0->0), fcn. (283->6), ass. (0->35)
t145 = -pkin(7) - pkin(6);
t132 = qJD(1) ^ 2;
t144 = t132 / 0.2e1;
t143 = cos(qJ(4));
t131 = cos(qJ(2));
t142 = t131 * t132;
t128 = sin(qJ(3));
t130 = cos(qJ(3));
t138 = qJD(1) * t131;
t129 = sin(qJ(2));
t139 = qJD(1) * t129;
t117 = t128 * t139 - t130 * t138;
t118 = (t128 * t131 + t129 * t130) * qJD(1);
t123 = (-pkin(2) * t131 - pkin(1)) * qJD(1);
t109 = t117 * pkin(3) - t118 * pkin(8) + t123;
t126 = qJD(2) + qJD(3);
t121 = qJD(2) * pkin(2) + t145 * t139;
t122 = t145 * t138;
t140 = t128 * t121 - t130 * t122;
t112 = t126 * pkin(8) + t140;
t127 = sin(qJ(4));
t141 = t127 * t109 + t143 * t112;
t137 = qJD(1) * qJD(2);
t136 = t129 * t137;
t135 = t131 * t137;
t134 = t143 * t109 - t127 * t112;
t133 = t130 * t121 + t128 * t122;
t111 = -t126 * pkin(3) - t133;
t116 = qJD(4) + t117;
t114 = t143 * t118 + t127 * t126;
t113 = t127 * t118 - t143 * t126;
t106 = t113 * pkin(4) + qJD(5) + t111;
t105 = -t113 * qJ(5) + t141;
t104 = t116 * pkin(4) - t114 * qJ(5) + t134;
t1 = [t144, 0, 0, t129 ^ 2 * t144, t129 * t142, t136, t135, qJD(2) ^ 2 / 0.2e1, pkin(1) * t142 - pkin(6) * t136, -t132 * pkin(1) * t129 - pkin(6) * t135, t118 ^ 2 / 0.2e1, -t118 * t117, t118 * t126, -t117 * t126, t126 ^ 2 / 0.2e1, t123 * t117 + t133 * t126, t123 * t118 - t140 * t126, t114 ^ 2 / 0.2e1, -t114 * t113, t114 * t116, -t113 * t116, t116 ^ 2 / 0.2e1, t111 * t113 + t134 * t116, t111 * t114 - t141 * t116, -t104 * t114 - t105 * t113, t105 ^ 2 / 0.2e1 + t104 ^ 2 / 0.2e1 + t106 ^ 2 / 0.2e1;];
T_reg = t1;
