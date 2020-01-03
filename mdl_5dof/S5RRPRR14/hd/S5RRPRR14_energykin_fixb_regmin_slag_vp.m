% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:36
% EndTime: 2019-12-31 20:38:37
% DurationCPUTime: 0.11s
% Computational Cost: add. (375->45), mult. (948->97), div. (0->0), fcn. (731->10), ass. (0->42)
t202 = cos(pkin(10));
t182 = sin(pkin(5));
t190 = qJD(1) ^ 2;
t201 = t182 ^ 2 * t190;
t197 = cos(pkin(5)) * qJD(1);
t179 = qJD(2) + t197;
t186 = sin(qJ(2));
t189 = cos(qJ(2));
t198 = qJD(1) * t182;
t193 = t189 * t198;
t196 = pkin(1) * t197;
t199 = pkin(7) * t193 + t186 * t196;
t167 = t179 * qJ(3) + t199;
t169 = (-pkin(2) * t189 - qJ(3) * t186 - pkin(1)) * t198;
t181 = sin(pkin(10));
t157 = -t181 * t167 + t202 * t169;
t194 = t186 * t198;
t171 = t181 * t179 + t202 * t194;
t152 = -pkin(3) * t193 - t171 * pkin(8) + t157;
t158 = t202 * t167 + t181 * t169;
t170 = -t202 * t179 + t181 * t194;
t154 = -t170 * pkin(8) + t158;
t185 = sin(qJ(4));
t188 = cos(qJ(4));
t200 = t185 * t152 + t188 * t154;
t195 = t189 * t201;
t161 = t188 * t170 + t185 * t171;
t192 = t188 * t152 - t185 * t154;
t191 = -pkin(7) * t194 + t189 * t196;
t164 = -t179 * pkin(2) + qJD(3) - t191;
t160 = t170 * pkin(3) + t164;
t187 = cos(qJ(5));
t184 = sin(qJ(5));
t174 = -qJD(4) + t193;
t162 = -t185 * t170 + t188 * t171;
t159 = qJD(5) + t161;
t156 = t187 * t162 - t184 * t174;
t155 = t184 * t162 + t187 * t174;
t150 = t161 * pkin(4) - t162 * pkin(9) + t160;
t149 = -t174 * pkin(9) + t200;
t148 = t174 * pkin(4) - t192;
t1 = [t190 / 0.2e1, 0, 0, t186 ^ 2 * t201 / 0.2e1, t186 * t195, t179 * t194, t179 * t193, t179 ^ 2 / 0.2e1, pkin(1) * t195 + t191 * t179, -pkin(1) * t186 * t201 - t199 * t179, -t157 * t193 + t164 * t170, t158 * t193 + t164 * t171, -t157 * t171 - t158 * t170, t158 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1, t162 ^ 2 / 0.2e1, -t162 * t161, -t162 * t174, t161 * t174, t174 ^ 2 / 0.2e1, t160 * t161 - t192 * t174, t160 * t162 + t200 * t174, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t159, -t155 * t159, t159 ^ 2 / 0.2e1, (-t184 * t149 + t187 * t150) * t159 + t148 * t155, -(t187 * t149 + t184 * t150) * t159 + t148 * t156;];
T_reg = t1;
