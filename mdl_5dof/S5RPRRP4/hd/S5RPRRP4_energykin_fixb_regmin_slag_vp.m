% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:06:51
% EndTime: 2019-12-05 18:06:51
% DurationCPUTime: 0.12s
% Computational Cost: add. (146->37), mult. (385->85), div. (0->0), fcn. (236->6), ass. (0->32)
t133 = sin(pkin(8));
t131 = t133 ^ 2;
t139 = qJD(1) ^ 2;
t149 = t131 * t139;
t134 = cos(pkin(8));
t121 = qJD(2) + (-pkin(2) * t134 - pkin(6) * t133 - pkin(1)) * qJD(1);
t138 = cos(qJ(3));
t120 = t138 * t121;
t145 = t134 * qJD(1);
t127 = -qJD(3) + t145;
t136 = sin(qJ(3));
t113 = -t127 * pkin(3) + t120 + (-pkin(7) * t133 * t138 - qJ(2) * t134 * t136) * qJD(1);
t146 = qJD(1) * t133;
t142 = t136 * t146;
t143 = qJ(2) * t145;
t147 = t136 * t121 + t138 * t143;
t116 = -pkin(7) * t142 + t147;
t135 = sin(qJ(4));
t137 = cos(qJ(4));
t148 = t135 * t113 + t137 * t116;
t122 = pkin(3) * t142 + qJ(2) * t146;
t144 = t136 * t149;
t141 = t137 * t113 - t135 * t116;
t132 = t134 ^ 2;
t130 = -qJD(1) * pkin(1) + qJD(2);
t124 = -qJD(4) + t127;
t118 = (-t135 * t136 + t137 * t138) * t146;
t117 = (t135 * t138 + t136 * t137) * t146;
t115 = t117 * pkin(4) + qJD(5) + t122;
t110 = -t117 * qJ(5) + t148;
t109 = -t124 * pkin(4) - t118 * qJ(5) + t141;
t1 = [t139 / 0.2e1, 0, 0, -t130 * t145, t130 * t146, (t131 + t132) * t139 * qJ(2), t130 ^ 2 / 0.2e1 + (t132 / 0.2e1 + t131 / 0.2e1) * qJ(2) ^ 2 * t139, t138 ^ 2 * t149 / 0.2e1, -t138 * t144, -t138 * t127 * t146, t127 * t142, t127 ^ 2 / 0.2e1, qJ(2) * t144 - (-t136 * t143 + t120) * t127, qJ(2) * t138 * t149 + t147 * t127, t118 ^ 2 / 0.2e1, -t118 * t117, -t118 * t124, t117 * t124, t124 ^ 2 / 0.2e1, t122 * t117 - t141 * t124, t122 * t118 + t148 * t124, -t109 * t118 - t110 * t117, t110 ^ 2 / 0.2e1 + t109 ^ 2 / 0.2e1 + t115 ^ 2 / 0.2e1;];
T_reg = t1;
