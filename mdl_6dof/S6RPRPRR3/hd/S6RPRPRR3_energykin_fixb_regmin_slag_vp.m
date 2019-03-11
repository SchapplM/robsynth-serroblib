% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR3
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
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:26
% EndTime: 2019-03-09 03:42:26
% DurationCPUTime: 0.12s
% Computational Cost: add. (336->48), mult. (753->105), div. (0->0), fcn. (508->10), ass. (0->42)
t185 = qJD(1) ^ 2;
t196 = t185 / 0.2e1;
t195 = cos(qJ(5));
t177 = sin(pkin(10));
t168 = (pkin(1) * t177 + pkin(7)) * qJD(1);
t182 = sin(qJ(3));
t184 = cos(qJ(3));
t193 = qJD(2) * t182 + t168 * t184;
t161 = qJD(3) * qJ(4) + t193;
t179 = cos(pkin(10));
t189 = -pkin(1) * t179 - pkin(2);
t162 = (-pkin(3) * t184 - qJ(4) * t182 + t189) * qJD(1);
t176 = sin(pkin(11));
t178 = cos(pkin(11));
t151 = -t161 * t176 + t162 * t178;
t192 = qJD(1) * t182;
t165 = qJD(3) * t176 + t178 * t192;
t191 = t184 * qJD(1);
t148 = -pkin(4) * t191 - pkin(8) * t165 + t151;
t152 = t161 * t178 + t162 * t176;
t164 = -qJD(3) * t178 + t176 * t192;
t150 = -pkin(8) * t164 + t152;
t181 = sin(qJ(5));
t194 = t148 * t181 + t150 * t195;
t190 = qJD(1) * qJD(3);
t188 = t148 * t195 - t150 * t181;
t187 = qJD(2) * t184 - t168 * t182;
t172 = -qJD(5) + t191;
t160 = -qJD(3) * pkin(3) + qJD(4) - t187;
t153 = pkin(4) * t164 + t160;
t183 = cos(qJ(6));
t180 = sin(qJ(6));
t170 = -qJD(6) + t172;
t169 = t189 * qJD(1);
t156 = -t164 * t181 + t165 * t195;
t155 = t164 * t195 + t165 * t181;
t147 = pkin(5) * t155 + t153;
t144 = -t155 * t180 + t156 * t183;
t143 = t155 * t183 + t156 * t180;
t142 = -pkin(9) * t155 + t194;
t141 = -pkin(5) * t172 - pkin(9) * t156 + t188;
t1 = [t196, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t177 ^ 2 / 0.2e1 + t179 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t185, t182 ^ 2 * t196, t182 * t185 * t184, t182 * t190, t184 * t190, qJD(3) ^ 2 / 0.2e1, qJD(3) * t187 - t169 * t191, -qJD(3) * t193 + t169 * t192, -t151 * t191 + t160 * t164, t152 * t191 + t160 * t165, -t151 * t165 - t152 * t164, t152 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1, t156 ^ 2 / 0.2e1, -t156 * t155, -t156 * t172, t155 * t172, t172 ^ 2 / 0.2e1, t153 * t155 - t172 * t188, t153 * t156 + t172 * t194, t144 ^ 2 / 0.2e1, -t144 * t143, -t144 * t170, t143 * t170, t170 ^ 2 / 0.2e1 -(t141 * t183 - t142 * t180) * t170 + t147 * t143 (t141 * t180 + t142 * t183) * t170 + t147 * t144;];
T_reg  = t1;
