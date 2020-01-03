% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:34
% EndTime: 2019-12-31 21:40:34
% DurationCPUTime: 0.12s
% Computational Cost: add. (416->45), mult. (1003->97), div. (0->0), fcn. (764->10), ass. (0->42)
t202 = cos(pkin(10));
t182 = sin(pkin(5));
t190 = qJD(1) ^ 2;
t201 = t182 ^ 2 * t190;
t197 = cos(pkin(5)) * qJD(1);
t179 = qJD(2) + t197;
t189 = cos(qJ(2));
t186 = sin(qJ(2));
t198 = qJD(1) * t182;
t194 = t186 * t198;
t196 = pkin(1) * t197;
t191 = -pkin(7) * t194 + t189 * t196;
t166 = -t179 * pkin(2) - t191;
t185 = sin(qJ(3));
t188 = cos(qJ(3));
t170 = -t188 * t179 + t185 * t194;
t171 = t185 * t179 + t188 * t194;
t156 = t170 * pkin(3) - t171 * qJ(4) + t166;
t193 = t189 * t198;
t174 = -qJD(3) + t193;
t199 = pkin(7) * t193 + t186 * t196;
t167 = t179 * pkin(8) + t199;
t168 = (-pkin(2) * t189 - pkin(8) * t186 - pkin(1)) * t198;
t200 = t188 * t167 + t185 * t168;
t159 = -t174 * qJ(4) + t200;
t181 = sin(pkin(10));
t150 = t181 * t156 + t202 * t159;
t195 = t189 * t201;
t149 = t202 * t156 - t181 * t159;
t192 = -t185 * t167 + t188 * t168;
t158 = t174 * pkin(3) + qJD(4) - t192;
t187 = cos(qJ(5));
t184 = sin(qJ(5));
t169 = qJD(5) + t170;
t162 = t202 * t171 - t181 * t174;
t161 = t181 * t171 + t202 * t174;
t153 = -t184 * t161 + t187 * t162;
t152 = t187 * t161 + t184 * t162;
t151 = t161 * pkin(4) + t158;
t148 = -t161 * pkin(9) + t150;
t147 = t170 * pkin(4) - t162 * pkin(9) + t149;
t1 = [t190 / 0.2e1, 0, 0, t186 ^ 2 * t201 / 0.2e1, t186 * t195, t179 * t194, t179 * t193, t179 ^ 2 / 0.2e1, pkin(1) * t195 + t191 * t179, -pkin(1) * t186 * t201 - t199 * t179, t171 ^ 2 / 0.2e1, -t171 * t170, -t171 * t174, t170 * t174, t174 ^ 2 / 0.2e1, t166 * t170 - t192 * t174, t166 * t171 + t200 * t174, t149 * t170 + t158 * t161, -t150 * t170 + t158 * t162, -t149 * t162 - t150 * t161, t150 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t169, -t152 * t169, t169 ^ 2 / 0.2e1, (t187 * t147 - t184 * t148) * t169 + t151 * t152, t151 * t153 - (t184 * t147 + t187 * t148) * t169;];
T_reg = t1;
