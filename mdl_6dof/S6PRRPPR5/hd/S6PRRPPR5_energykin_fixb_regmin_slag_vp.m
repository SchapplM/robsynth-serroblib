% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:05
% EndTime: 2019-03-08 21:21:05
% DurationCPUTime: 0.10s
% Computational Cost: add. (232->47), mult. (511->98), div. (0->0), fcn. (330->10), ass. (0->41)
t186 = qJD(2) ^ 2;
t199 = t186 / 0.2e1;
t198 = -pkin(3) - qJ(5);
t182 = sin(qJ(2));
t196 = qJD(1) * sin(pkin(6));
t170 = qJD(2) * pkin(8) + t182 * t196;
t181 = sin(qJ(3));
t184 = cos(qJ(3));
t195 = qJD(1) * cos(pkin(6));
t188 = -t181 * t170 + t184 * t195;
t187 = qJD(4) - t188;
t193 = t181 * qJD(2);
t156 = pkin(4) * t193 + t198 * qJD(3) + t187;
t189 = -qJ(4) * t181 - pkin(2);
t185 = cos(qJ(2));
t191 = t185 * t196;
t161 = -t191 + (t198 * t184 + t189) * qJD(2);
t176 = sin(pkin(11));
t178 = cos(pkin(11));
t152 = t176 * t156 + t178 * t161;
t197 = t184 * t170 + t181 * t195;
t194 = qJD(2) * t184;
t192 = qJD(2) * qJD(3);
t163 = -qJD(3) * qJ(4) - t197;
t190 = qJD(2) * t196;
t151 = t178 * t156 - t176 * t161;
t159 = pkin(4) * t194 + qJD(5) - t163;
t183 = cos(qJ(6));
t180 = sin(qJ(6));
t173 = qJD(6) + t193;
t171 = -qJD(2) * pkin(2) - t191;
t169 = t178 * qJD(3) - t176 * t194;
t168 = t176 * qJD(3) + t178 * t194;
t164 = -t191 + (-pkin(3) * t184 + t189) * qJD(2);
t162 = -qJD(3) * pkin(3) + t187;
t158 = -t180 * t168 + t183 * t169;
t157 = t183 * t168 + t180 * t169;
t153 = t168 * pkin(5) + t159;
t150 = -t168 * pkin(9) + t152;
t149 = pkin(5) * t193 - t169 * pkin(9) + t151;
t1 = [qJD(1) ^ 2 / 0.2e1, t199, t185 * t190, -t182 * t190, t181 ^ 2 * t199, t181 * t186 * t184, t181 * t192, t184 * t192, qJD(3) ^ 2 / 0.2e1, t188 * qJD(3) - t171 * t194, -t197 * qJD(3) + t171 * t193 (t162 * t181 - t163 * t184) * qJD(2), t162 * qJD(3) + t164 * t194, -t163 * qJD(3) - t164 * t193, t164 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, t151 * t193 + t159 * t168, -t152 * t193 + t159 * t169, -t151 * t169 - t152 * t168, t152 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1, t158 ^ 2 / 0.2e1, -t158 * t157, t158 * t173, -t157 * t173, t173 ^ 2 / 0.2e1 (t183 * t149 - t180 * t150) * t173 + t153 * t157 -(t180 * t149 + t183 * t150) * t173 + t153 * t158;];
T_reg  = t1;
