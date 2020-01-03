% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP7
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
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:43
% EndTime: 2019-12-31 21:57:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (251->37), mult. (559->82), div. (0->0), fcn. (371->6), ass. (0->35)
t153 = -pkin(7) - pkin(6);
t141 = qJD(1) ^ 2;
t152 = t141 / 0.2e1;
t140 = cos(qJ(2));
t151 = t140 * t141;
t136 = sin(qJ(3));
t139 = cos(qJ(3));
t147 = qJD(1) * t140;
t137 = sin(qJ(2));
t148 = qJD(1) * t137;
t125 = t136 * t148 - t139 * t147;
t126 = (t136 * t140 + t137 * t139) * qJD(1);
t131 = (-pkin(2) * t140 - pkin(1)) * qJD(1);
t117 = t125 * pkin(3) - t126 * pkin(8) + t131;
t134 = qJD(2) + qJD(3);
t129 = qJD(2) * pkin(2) + t153 * t148;
t130 = t153 * t147;
t149 = t136 * t129 - t139 * t130;
t120 = t134 * pkin(8) + t149;
t135 = sin(qJ(4));
t138 = cos(qJ(4));
t150 = t135 * t117 + t138 * t120;
t146 = qJD(1) * qJD(2);
t145 = t137 * t146;
t144 = t140 * t146;
t143 = t139 * t129 + t136 * t130;
t142 = t138 * t117 - t135 * t120;
t119 = -t134 * pkin(3) - t143;
t124 = qJD(4) + t125;
t122 = t138 * t126 + t135 * t134;
t121 = t135 * t126 - t138 * t134;
t115 = t121 * pkin(4) - t122 * qJ(5) + t119;
t114 = t124 * qJ(5) + t150;
t113 = -t124 * pkin(4) + qJD(5) - t142;
t1 = [t152, 0, 0, t137 ^ 2 * t152, t137 * t151, t145, t144, qJD(2) ^ 2 / 0.2e1, pkin(1) * t151 - pkin(6) * t145, -t141 * pkin(1) * t137 - pkin(6) * t144, t126 ^ 2 / 0.2e1, -t126 * t125, t126 * t134, -t125 * t134, t134 ^ 2 / 0.2e1, t131 * t125 + t143 * t134, t131 * t126 - t149 * t134, t122 ^ 2 / 0.2e1, -t122 * t121, t122 * t124, -t121 * t124, t124 ^ 2 / 0.2e1, t119 * t121 + t142 * t124, t119 * t122 - t150 * t124, -t113 * t124 + t115 * t121, t113 * t122 - t114 * t121, t114 * t124 - t115 * t122, t114 ^ 2 / 0.2e1 + t115 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1;];
T_reg = t1;
