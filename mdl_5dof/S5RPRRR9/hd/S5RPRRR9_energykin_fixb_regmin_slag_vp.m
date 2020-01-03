% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:59
% EndTime: 2019-12-31 19:07:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (225->42), mult. (604->88), div. (0->0), fcn. (448->8), ass. (0->37)
t156 = cos(qJ(3));
t155 = pkin(6) + qJ(2);
t140 = sin(pkin(9));
t141 = cos(pkin(9));
t144 = sin(qJ(3));
t130 = (t156 * t140 + t141 * t144) * qJD(1);
t152 = qJD(1) * t140;
t131 = t155 * t152;
t151 = qJD(1) * t141;
t132 = t155 * t151;
t150 = -t156 * t131 - t132 * t144;
t117 = qJD(3) * pkin(3) - pkin(7) * t130 + t150;
t129 = t144 * t152 - t156 * t151;
t153 = -t144 * t131 + t156 * t132;
t118 = -pkin(7) * t129 + t153;
t143 = sin(qJ(4));
t146 = cos(qJ(4));
t154 = t143 * t117 + t146 * t118;
t122 = t146 * t129 + t130 * t143;
t149 = t117 * t146 - t118 * t143;
t133 = qJD(2) + (-pkin(2) * t141 - pkin(1)) * qJD(1);
t124 = pkin(3) * t129 + t133;
t147 = qJD(1) ^ 2;
t145 = cos(qJ(5));
t142 = sin(qJ(5));
t139 = qJD(3) + qJD(4);
t138 = t141 ^ 2;
t137 = t140 ^ 2;
t136 = -qJD(1) * pkin(1) + qJD(2);
t123 = -t129 * t143 + t130 * t146;
t121 = qJD(5) + t122;
t120 = t123 * t145 + t139 * t142;
t119 = t123 * t142 - t145 * t139;
t114 = pkin(4) * t122 - pkin(8) * t123 + t124;
t113 = pkin(8) * t139 + t154;
t112 = -pkin(4) * t139 - t149;
t1 = [t147 / 0.2e1, 0, 0, -t136 * t151, t136 * t152, (t137 + t138) * t147 * qJ(2), t136 ^ 2 / 0.2e1 + (t138 / 0.2e1 + t137 / 0.2e1) * qJ(2) ^ 2 * t147, t130 ^ 2 / 0.2e1, -t130 * t129, t130 * qJD(3), -t129 * qJD(3), qJD(3) ^ 2 / 0.2e1, t150 * qJD(3) + t133 * t129, -t153 * qJD(3) + t133 * t130, t123 ^ 2 / 0.2e1, -t123 * t122, t123 * t139, -t122 * t139, t139 ^ 2 / 0.2e1, t124 * t122 + t149 * t139, t124 * t123 - t154 * t139, t120 ^ 2 / 0.2e1, -t120 * t119, t120 * t121, -t119 * t121, t121 ^ 2 / 0.2e1, (-t142 * t113 + t114 * t145) * t121 + t112 * t119, -(t113 * t145 + t142 * t114) * t121 + t112 * t120;];
T_reg = t1;
