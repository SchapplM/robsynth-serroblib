% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:38:57
% EndTime: 2019-03-09 03:38:57
% DurationCPUTime: 0.12s
% Computational Cost: add. (292->47), mult. (680->103), div. (0->0), fcn. (467->10), ass. (0->42)
t188 = qJD(1) ^ 2;
t198 = t188 / 0.2e1;
t197 = cos(qJ(5));
t180 = sin(pkin(10));
t173 = (pkin(1) * t180 + pkin(7)) * qJD(1);
t187 = cos(qJ(3));
t178 = t187 * qJD(2);
t185 = sin(qJ(3));
t163 = qJD(3) * pkin(3) + t178 + (-qJ(4) * qJD(1) - t173) * t185;
t193 = qJD(1) * t187;
t195 = qJD(2) * t185 + t173 * t187;
t166 = qJ(4) * t193 + t195;
t179 = sin(pkin(11));
t181 = cos(pkin(11));
t153 = t163 * t179 + t166 * t181;
t151 = qJD(3) * pkin(8) + t153;
t182 = cos(pkin(10));
t191 = -pkin(1) * t182 - pkin(2);
t169 = qJD(4) + (-pkin(3) * t187 + t191) * qJD(1);
t194 = qJD(1) * t185;
t170 = -t179 * t194 + t181 * t193;
t171 = (t179 * t187 + t181 * t185) * qJD(1);
t158 = -t170 * pkin(4) - t171 * pkin(8) + t169;
t184 = sin(qJ(5));
t196 = t151 * t197 + t158 * t184;
t192 = qJD(1) * qJD(3);
t190 = -t151 * t184 + t158 * t197;
t152 = t163 * t181 - t166 * t179;
t168 = qJD(5) - t170;
t150 = -qJD(3) * pkin(4) - t152;
t186 = cos(qJ(6));
t183 = sin(qJ(6));
t174 = t191 * qJD(1);
t167 = qJD(6) + t168;
t165 = qJD(3) * t184 + t171 * t197;
t164 = -qJD(3) * t197 + t171 * t184;
t155 = -t164 * t183 + t165 * t186;
t154 = t164 * t186 + t165 * t183;
t148 = pkin(5) * t164 + t150;
t147 = -pkin(9) * t164 + t196;
t146 = pkin(5) * t168 - pkin(9) * t165 + t190;
t1 = [t198, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t180 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t188, t185 ^ 2 * t198, t185 * t188 * t187, t185 * t192, t187 * t192, qJD(3) ^ 2 / 0.2e1, -t174 * t193 + (-t173 * t185 + t178) * qJD(3), -qJD(3) * t195 + t174 * t194, -t152 * t171 + t153 * t170, t153 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1, t165 ^ 2 / 0.2e1, -t165 * t164, t165 * t168, -t164 * t168, t168 ^ 2 / 0.2e1, t150 * t164 + t168 * t190, t150 * t165 - t168 * t196, t155 ^ 2 / 0.2e1, -t155 * t154, t155 * t167, -t154 * t167, t167 ^ 2 / 0.2e1 (t146 * t186 - t147 * t183) * t167 + t148 * t154 -(t146 * t183 + t147 * t186) * t167 + t148 * t155;];
T_reg  = t1;
