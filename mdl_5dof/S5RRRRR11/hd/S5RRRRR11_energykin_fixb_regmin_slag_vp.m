% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:41
% EndTime: 2019-12-31 22:43:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (358->44), mult. (851->96), div. (0->0), fcn. (668->10), ass. (0->43)
t204 = cos(qJ(4));
t181 = sin(pkin(5));
t190 = qJD(1) ^ 2;
t203 = t181 ^ 2 * t190;
t198 = cos(pkin(5)) * qJD(1);
t179 = qJD(2) + t198;
t189 = cos(qJ(2));
t186 = sin(qJ(2));
t199 = qJD(1) * t181;
t195 = t186 * t199;
t197 = pkin(1) * t198;
t191 = -pkin(7) * t195 + t189 * t197;
t165 = -t179 * pkin(2) - t191;
t185 = sin(qJ(3));
t188 = cos(qJ(3));
t170 = -t188 * t179 + t185 * t195;
t171 = t185 * t179 + t188 * t195;
t155 = t170 * pkin(3) - t171 * pkin(9) + t165;
t194 = t189 * t199;
t174 = -qJD(3) + t194;
t200 = pkin(7) * t194 + t186 * t197;
t166 = t179 * pkin(8) + t200;
t168 = (-pkin(2) * t189 - pkin(8) * t186 - pkin(1)) * t199;
t201 = t188 * t166 + t185 * t168;
t158 = -t174 * pkin(9) + t201;
t184 = sin(qJ(4));
t202 = t184 * t155 + t204 * t158;
t196 = t189 * t203;
t193 = t204 * t155 - t184 * t158;
t192 = -t185 * t166 + t188 * t168;
t157 = t174 * pkin(3) - t192;
t169 = qJD(4) + t170;
t187 = cos(qJ(5));
t183 = sin(qJ(5));
t167 = qJD(5) + t169;
t161 = t204 * t171 - t184 * t174;
t160 = t184 * t171 + t204 * t174;
t152 = -t183 * t160 + t187 * t161;
t151 = t187 * t160 + t183 * t161;
t150 = t160 * pkin(4) + t157;
t149 = -t160 * pkin(10) + t202;
t148 = t169 * pkin(4) - t161 * pkin(10) + t193;
t1 = [t190 / 0.2e1, 0, 0, t186 ^ 2 * t203 / 0.2e1, t186 * t196, t179 * t195, t179 * t194, t179 ^ 2 / 0.2e1, pkin(1) * t196 + t191 * t179, -pkin(1) * t186 * t203 - t200 * t179, t171 ^ 2 / 0.2e1, -t171 * t170, -t171 * t174, t170 * t174, t174 ^ 2 / 0.2e1, t165 * t170 - t192 * t174, t165 * t171 + t201 * t174, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t169, -t160 * t169, t169 ^ 2 / 0.2e1, t157 * t160 + t193 * t169, t157 * t161 - t202 * t169, t152 ^ 2 / 0.2e1, -t152 * t151, t152 * t167, -t151 * t167, t167 ^ 2 / 0.2e1, t150 * t151 + (t187 * t148 - t183 * t149) * t167, t150 * t152 - (t183 * t148 + t187 * t149) * t167;];
T_reg = t1;
