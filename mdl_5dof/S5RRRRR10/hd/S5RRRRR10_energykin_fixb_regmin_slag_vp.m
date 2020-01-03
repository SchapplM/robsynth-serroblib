% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR10
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:35:57
% EndTime: 2019-12-31 22:35:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (371->44), mult. (888->96), div. (0->0), fcn. (699->10), ass. (0->43)
t206 = cos(qJ(3));
t183 = sin(pkin(5));
t192 = qJD(1) ^ 2;
t205 = t183 ^ 2 * t192;
t200 = cos(pkin(5)) * qJD(1);
t181 = qJD(2) + t200;
t187 = sin(qJ(3));
t188 = sin(qJ(2));
t201 = qJD(1) * t183;
t197 = t188 * t201;
t171 = t187 * t181 + t206 * t197;
t191 = cos(qJ(2));
t196 = t191 * t201;
t176 = -qJD(3) + t196;
t199 = pkin(1) * t200;
t202 = pkin(7) * t196 + t188 * t199;
t167 = t181 * pkin(8) + t202;
t169 = (-pkin(2) * t191 - pkin(8) * t188 - pkin(1)) * t201;
t195 = -t187 * t167 + t206 * t169;
t154 = -t176 * pkin(3) - t171 * pkin(9) + t195;
t170 = -t206 * t181 + t187 * t197;
t203 = t206 * t167 + t187 * t169;
t156 = -t170 * pkin(9) + t203;
t186 = sin(qJ(4));
t190 = cos(qJ(4));
t204 = t186 * t154 + t190 * t156;
t198 = t191 * t205;
t160 = t190 * t170 + t186 * t171;
t194 = t190 * t154 - t186 * t156;
t193 = -pkin(7) * t197 + t191 * t199;
t166 = -t181 * pkin(2) - t193;
t162 = t170 * pkin(3) + t166;
t189 = cos(qJ(5));
t185 = sin(qJ(5));
t173 = -qJD(4) + t176;
t161 = -t186 * t170 + t190 * t171;
t159 = qJD(5) + t160;
t158 = t189 * t161 - t185 * t173;
t157 = t185 * t161 + t189 * t173;
t152 = t160 * pkin(4) - t161 * pkin(10) + t162;
t151 = -t173 * pkin(10) + t204;
t150 = t173 * pkin(4) - t194;
t1 = [t192 / 0.2e1, 0, 0, t188 ^ 2 * t205 / 0.2e1, t188 * t198, t181 * t197, t181 * t196, t181 ^ 2 / 0.2e1, pkin(1) * t198 + t193 * t181, -pkin(1) * t188 * t205 - t202 * t181, t171 ^ 2 / 0.2e1, -t171 * t170, -t171 * t176, t170 * t176, t176 ^ 2 / 0.2e1, t166 * t170 - t195 * t176, t166 * t171 + t203 * t176, t161 ^ 2 / 0.2e1, -t161 * t160, -t161 * t173, t160 * t173, t173 ^ 2 / 0.2e1, t162 * t160 - t194 * t173, t162 * t161 + t204 * t173, t158 ^ 2 / 0.2e1, -t158 * t157, t158 * t159, -t157 * t159, t159 ^ 2 / 0.2e1, (-t185 * t151 + t189 * t152) * t159 + t150 * t157, -(t189 * t151 + t185 * t152) * t159 + t150 * t158;];
T_reg = t1;
