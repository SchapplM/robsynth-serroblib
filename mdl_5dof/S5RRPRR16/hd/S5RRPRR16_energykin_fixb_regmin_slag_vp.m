% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR16_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:25
% EndTime: 2019-12-31 20:47:25
% DurationCPUTime: 0.12s
% Computational Cost: add. (213->44), mult. (539->95), div. (0->0), fcn. (376->8), ass. (0->40)
t179 = -pkin(2) - pkin(8);
t158 = sin(pkin(5));
t166 = qJD(1) ^ 2;
t178 = t158 ^ 2 * t166;
t162 = sin(qJ(2));
t175 = qJD(1) * t158;
t170 = t162 * t175;
t152 = pkin(7) * t170;
t159 = cos(pkin(5));
t174 = t159 * qJD(1);
t156 = qJD(2) + t174;
t165 = cos(qJ(2));
t137 = qJD(3) + t152 + t179 * t156 + (-pkin(1) * t159 * t165 + pkin(3) * t158 * t162) * qJD(1);
t169 = -qJ(3) * t162 - pkin(1);
t142 = (t179 * t165 + t169) * t175;
t161 = sin(qJ(4));
t164 = cos(qJ(4));
t177 = t161 * t137 + t164 * t142;
t171 = t165 * t175;
t173 = pkin(1) * t174;
t176 = pkin(7) * t171 + t162 * t173;
t172 = t165 * t178;
t144 = -t156 * qJ(3) - t176;
t141 = pkin(3) * t171 - t144;
t168 = t164 * t137 - t161 * t142;
t167 = t165 * t173 - t152;
t147 = t161 * t156 + t164 * t171;
t163 = cos(qJ(5));
t160 = sin(qJ(5));
t150 = qJD(4) + t170;
t148 = t164 * t156 - t161 * t171;
t146 = qJD(5) + t147;
t145 = (-pkin(2) * t165 + t169) * t175;
t143 = -t156 * pkin(2) + qJD(3) - t167;
t139 = t163 * t148 + t160 * t150;
t138 = t160 * t148 - t163 * t150;
t135 = t147 * pkin(4) - t148 * pkin(9) + t141;
t134 = t150 * pkin(9) + t177;
t133 = -t150 * pkin(4) - t168;
t1 = [t166 / 0.2e1, 0, 0, t162 ^ 2 * t178 / 0.2e1, t162 * t172, t156 * t170, t156 * t171, t156 ^ 2 / 0.2e1, pkin(1) * t172 + t167 * t156, -pkin(1) * t162 * t178 - t176 * t156, (t143 * t162 - t144 * t165) * t175, t143 * t156 + t145 * t171, -t144 * t156 - t145 * t170, t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1, t148 ^ 2 / 0.2e1, -t148 * t147, t148 * t150, -t147 * t150, t150 ^ 2 / 0.2e1, t141 * t147 + t168 * t150, t141 * t148 - t177 * t150, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * t146, -t138 * t146, t146 ^ 2 / 0.2e1, (-t160 * t134 + t163 * t135) * t146 + t133 * t138, -(t163 * t134 + t160 * t135) * t146 + t133 * t139;];
T_reg = t1;
