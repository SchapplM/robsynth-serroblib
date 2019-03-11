% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:29:47
% EndTime: 2019-03-10 03:29:47
% DurationCPUTime: 0.13s
% Computational Cost: add. (610->53), mult. (1461->112), div. (0->0), fcn. (1155->10), ass. (0->50)
t202 = -pkin(8) - pkin(7);
t186 = qJD(1) ^ 2;
t201 = t186 / 0.2e1;
t200 = cos(qJ(3));
t199 = cos(qJ(4));
t185 = cos(qJ(2));
t198 = t185 * t186;
t181 = sin(qJ(3));
t193 = qJD(1) * t185;
t182 = sin(qJ(2));
t194 = qJD(1) * t182;
t167 = t181 * t194 - t200 * t193;
t168 = (t181 * t185 + t200 * t182) * qJD(1);
t180 = sin(qJ(4));
t162 = -t180 * t167 + t199 * t168;
t177 = qJD(2) + qJD(3);
t176 = qJD(4) + t177;
t170 = qJD(2) * pkin(2) + t202 * t194;
t171 = t202 * t193;
t188 = t200 * t170 + t181 * t171;
t157 = t177 * pkin(3) - t168 * pkin(9) + t188;
t195 = t181 * t170 - t200 * t171;
t159 = -t167 * pkin(9) + t195;
t189 = t199 * t157 - t180 * t159;
t146 = t176 * pkin(4) - t162 * pkin(10) + t189;
t161 = t199 * t167 + t180 * t168;
t196 = t180 * t157 + t199 * t159;
t148 = -t161 * pkin(10) + t196;
t179 = sin(qJ(5));
t184 = cos(qJ(5));
t197 = t179 * t146 + t184 * t148;
t192 = qJD(1) * qJD(2);
t191 = t182 * t192;
t190 = t185 * t192;
t152 = t184 * t161 + t179 * t162;
t173 = (-pkin(2) * t185 - pkin(1)) * qJD(1);
t187 = t184 * t146 - t179 * t148;
t163 = t167 * pkin(3) + t173;
t154 = t161 * pkin(4) + t163;
t183 = cos(qJ(6));
t178 = sin(qJ(6));
t175 = qJD(5) + t176;
t153 = -t179 * t161 + t184 * t162;
t151 = qJD(6) + t152;
t150 = t183 * t153 + t178 * t175;
t149 = t178 * t153 - t183 * t175;
t144 = t152 * pkin(5) - t153 * pkin(11) + t154;
t143 = t175 * pkin(11) + t197;
t142 = -t175 * pkin(5) - t187;
t1 = [t201, 0, 0, t182 ^ 2 * t201, t182 * t198, t191, t190, qJD(2) ^ 2 / 0.2e1, pkin(1) * t198 - pkin(7) * t191, -t186 * pkin(1) * t182 - pkin(7) * t190, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t177, -t167 * t177, t177 ^ 2 / 0.2e1, t173 * t167 + t188 * t177, t173 * t168 - t195 * t177, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * t176, -t161 * t176, t176 ^ 2 / 0.2e1, t163 * t161 + t189 * t176, t163 * t162 - t196 * t176, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t175, -t152 * t175, t175 ^ 2 / 0.2e1, t154 * t152 + t187 * t175, t154 * t153 - t197 * t175, t150 ^ 2 / 0.2e1, -t150 * t149, t150 * t151, -t149 * t151, t151 ^ 2 / 0.2e1 (-t178 * t143 + t183 * t144) * t151 + t142 * t149 -(t183 * t143 + t178 * t144) * t151 + t142 * t150;];
T_reg  = t1;
