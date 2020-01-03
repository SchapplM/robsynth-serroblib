% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:26
% EndTime: 2019-12-31 22:22:26
% DurationCPUTime: 0.13s
% Computational Cost: add. (255->40), mult. (614->89), div. (0->0), fcn. (451->8), ass. (0->40)
t159 = -pkin(7) - pkin(6);
t146 = qJD(1) ^ 2;
t158 = t146 / 0.2e1;
t157 = cos(qJ(3));
t145 = cos(qJ(2));
t156 = t145 * t146;
t141 = sin(qJ(3));
t142 = sin(qJ(2));
t130 = (t141 * t145 + t157 * t142) * qJD(1);
t138 = qJD(2) + qJD(3);
t153 = qJD(1) * t142;
t132 = qJD(2) * pkin(2) + t159 * t153;
t152 = qJD(1) * t145;
t133 = t159 * t152;
t148 = t157 * t132 + t141 * t133;
t117 = t138 * pkin(3) - t130 * pkin(8) + t148;
t129 = t141 * t153 - t157 * t152;
t154 = t141 * t132 - t157 * t133;
t119 = -t129 * pkin(8) + t154;
t140 = sin(qJ(4));
t144 = cos(qJ(4));
t155 = t140 * t117 + t144 * t119;
t151 = qJD(1) * qJD(2);
t150 = t142 * t151;
t149 = t145 * t151;
t123 = t144 * t129 + t140 * t130;
t134 = (-pkin(2) * t145 - pkin(1)) * qJD(1);
t147 = t144 * t117 - t140 * t119;
t125 = t129 * pkin(3) + t134;
t143 = cos(qJ(5));
t139 = sin(qJ(5));
t137 = qJD(4) + t138;
t124 = -t140 * t129 + t144 * t130;
t122 = qJD(5) + t123;
t121 = t143 * t124 + t139 * t137;
t120 = t139 * t124 - t143 * t137;
t115 = t123 * pkin(4) - t124 * pkin(9) + t125;
t114 = t137 * pkin(9) + t155;
t113 = -t137 * pkin(4) - t147;
t1 = [t158, 0, 0, t142 ^ 2 * t158, t142 * t156, t150, t149, qJD(2) ^ 2 / 0.2e1, pkin(1) * t156 - pkin(6) * t150, -t146 * pkin(1) * t142 - pkin(6) * t149, t130 ^ 2 / 0.2e1, -t130 * t129, t130 * t138, -t129 * t138, t138 ^ 2 / 0.2e1, t134 * t129 + t148 * t138, t134 * t130 - t154 * t138, t124 ^ 2 / 0.2e1, -t124 * t123, t124 * t137, -t123 * t137, t137 ^ 2 / 0.2e1, t125 * t123 + t147 * t137, t125 * t124 - t155 * t137, t121 ^ 2 / 0.2e1, -t121 * t120, t121 * t122, -t120 * t122, t122 ^ 2 / 0.2e1, (-t139 * t114 + t143 * t115) * t122 + t113 * t120, -(t143 * t114 + t139 * t115) * t122 + t113 * t121;];
T_reg = t1;
