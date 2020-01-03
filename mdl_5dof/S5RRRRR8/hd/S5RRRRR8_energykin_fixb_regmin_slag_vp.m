% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR8
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:46
% EndTime: 2019-12-31 22:25:46
% DurationCPUTime: 0.09s
% Computational Cost: add. (253->40), mult. (562->89), div. (0->0), fcn. (405->8), ass. (0->40)
t167 = -pkin(7) - pkin(6);
t154 = qJD(1) ^ 2;
t166 = t154 / 0.2e1;
t165 = cos(qJ(4));
t153 = cos(qJ(2));
t164 = t153 * t154;
t149 = sin(qJ(3));
t152 = cos(qJ(3));
t160 = qJD(1) * t153;
t150 = sin(qJ(2));
t161 = qJD(1) * t150;
t137 = t149 * t161 - t152 * t160;
t138 = (t149 * t153 + t150 * t152) * qJD(1);
t143 = (-pkin(2) * t153 - pkin(1)) * qJD(1);
t127 = t137 * pkin(3) - t138 * pkin(8) + t143;
t146 = qJD(2) + qJD(3);
t141 = qJD(2) * pkin(2) + t167 * t161;
t142 = t167 * t160;
t162 = t149 * t141 - t152 * t142;
t130 = t146 * pkin(8) + t162;
t148 = sin(qJ(4));
t163 = t148 * t127 + t165 * t130;
t159 = qJD(1) * qJD(2);
t158 = t150 * t159;
t157 = t153 * t159;
t156 = t165 * t127 - t148 * t130;
t155 = t152 * t141 + t149 * t142;
t129 = -t146 * pkin(3) - t155;
t136 = qJD(4) + t137;
t151 = cos(qJ(5));
t147 = sin(qJ(5));
t134 = qJD(5) + t136;
t133 = t165 * t138 + t148 * t146;
t132 = t148 * t138 - t165 * t146;
t124 = -t147 * t132 + t151 * t133;
t123 = t151 * t132 + t147 * t133;
t122 = t132 * pkin(4) + t129;
t121 = -t132 * pkin(9) + t163;
t120 = t136 * pkin(4) - t133 * pkin(9) + t156;
t1 = [t166, 0, 0, t150 ^ 2 * t166, t150 * t164, t158, t157, qJD(2) ^ 2 / 0.2e1, pkin(1) * t164 - pkin(6) * t158, -t154 * pkin(1) * t150 - pkin(6) * t157, t138 ^ 2 / 0.2e1, -t138 * t137, t138 * t146, -t137 * t146, t146 ^ 2 / 0.2e1, t143 * t137 + t155 * t146, t143 * t138 - t162 * t146, t133 ^ 2 / 0.2e1, -t133 * t132, t133 * t136, -t132 * t136, t136 ^ 2 / 0.2e1, t129 * t132 + t156 * t136, t129 * t133 - t163 * t136, t124 ^ 2 / 0.2e1, -t124 * t123, t124 * t134, -t123 * t134, t134 ^ 2 / 0.2e1, (t151 * t120 - t147 * t121) * t134 + t122 * t123, -(t147 * t120 + t151 * t121) * t134 + t122 * t124;];
T_reg = t1;
