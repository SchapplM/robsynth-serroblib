% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:04
% EndTime: 2019-12-31 19:33:04
% DurationCPUTime: 0.13s
% Computational Cost: add. (261->40), mult. (649->87), div. (0->0), fcn. (440->8), ass. (0->38)
t163 = qJD(1) ^ 2;
t172 = t163 / 0.2e1;
t171 = pkin(6) + qJ(3);
t170 = cos(pkin(9));
t162 = cos(qJ(2));
t169 = t162 * t163;
t157 = sin(pkin(8));
t158 = cos(pkin(8));
t167 = qJD(1) * t162;
t160 = sin(qJ(2));
t168 = qJD(1) * t160;
t147 = t157 * t168 - t158 * t167;
t148 = (t157 * t162 + t158 * t160) * qJD(1);
t153 = qJD(3) + (-pkin(2) * t162 - pkin(1)) * qJD(1);
t136 = t147 * pkin(3) - t148 * qJ(4) + t153;
t151 = qJD(2) * pkin(2) - t171 * t168;
t152 = t171 * t167;
t141 = t157 * t151 + t158 * t152;
t139 = qJD(2) * qJ(4) + t141;
t156 = sin(pkin(9));
t130 = t156 * t136 + t170 * t139;
t166 = qJD(1) * qJD(2);
t165 = t160 * t166;
t164 = t162 * t166;
t129 = t170 * t136 - t156 * t139;
t140 = t158 * t151 - t157 * t152;
t138 = -qJD(2) * pkin(3) + qJD(4) - t140;
t161 = cos(qJ(5));
t159 = sin(qJ(5));
t146 = qJD(5) + t147;
t144 = t156 * qJD(2) + t170 * t148;
t143 = -t170 * qJD(2) + t156 * t148;
t133 = -t159 * t143 + t161 * t144;
t132 = t161 * t143 + t159 * t144;
t131 = t143 * pkin(4) + t138;
t128 = -t143 * pkin(7) + t130;
t127 = t147 * pkin(4) - t144 * pkin(7) + t129;
t1 = [t172, 0, 0, t160 ^ 2 * t172, t160 * t169, t165, t164, qJD(2) ^ 2 / 0.2e1, pkin(1) * t169 - pkin(6) * t165, -t163 * pkin(1) * t160 - pkin(6) * t164, -t140 * t148 - t141 * t147, t141 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1, t129 * t147 + t138 * t143, -t130 * t147 + t138 * t144, -t129 * t144 - t130 * t143, t130 ^ 2 / 0.2e1 + t129 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1, t133 ^ 2 / 0.2e1, -t133 * t132, t133 * t146, -t132 * t146, t146 ^ 2 / 0.2e1, (t161 * t127 - t159 * t128) * t146 + t131 * t132, -(t159 * t127 + t161 * t128) * t146 + t131 * t133;];
T_reg = t1;
