% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR10
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
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:38
% EndTime: 2019-12-31 20:26:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (303->44), mult. (880->94), div. (0->0), fcn. (684->10), ass. (0->42)
t188 = sin(pkin(5));
t197 = qJD(1) ^ 2;
t207 = t188 ^ 2 * t197;
t196 = cos(qJ(2));
t203 = cos(pkin(5)) * qJD(1);
t202 = pkin(1) * t203;
t184 = t196 * t202;
t185 = qJD(2) + t203;
t193 = sin(qJ(2));
t204 = qJD(1) * t188;
t200 = t193 * t204;
t171 = t185 * pkin(2) + t184 + (-pkin(7) - qJ(3)) * t200;
t199 = t196 * t204;
t205 = pkin(7) * t199 + t193 * t202;
t174 = qJ(3) * t199 + t205;
t187 = sin(pkin(10));
t189 = cos(pkin(10));
t162 = t187 * t171 + t189 * t174;
t160 = t185 * pkin(8) + t162;
t177 = -t187 * t200 + t189 * t199;
t178 = (t187 * t196 + t189 * t193) * t204;
t179 = qJD(3) + (-pkin(2) * t196 - pkin(1)) * t204;
t166 = -t177 * pkin(3) - t178 * pkin(8) + t179;
t192 = sin(qJ(4));
t195 = cos(qJ(4));
t206 = t195 * t160 + t192 * t166;
t201 = t196 * t207;
t161 = t189 * t171 - t187 * t174;
t169 = t192 * t178 - t195 * t185;
t198 = -t192 * t160 + t195 * t166;
t159 = -t185 * pkin(3) - t161;
t194 = cos(qJ(5));
t191 = sin(qJ(5));
t176 = qJD(4) - t177;
t170 = t195 * t178 + t192 * t185;
t168 = qJD(5) + t169;
t164 = t194 * t170 + t191 * t176;
t163 = t191 * t170 - t194 * t176;
t157 = t169 * pkin(4) - t170 * pkin(9) + t159;
t156 = t176 * pkin(9) + t206;
t155 = -t176 * pkin(4) - t198;
t1 = [t197 / 0.2e1, 0, 0, t193 ^ 2 * t207 / 0.2e1, t193 * t201, t185 * t200, t185 * t199, t185 ^ 2 / 0.2e1, pkin(1) * t201 + (-pkin(7) * t200 + t184) * t185, -pkin(1) * t193 * t207 - t205 * t185, -t161 * t178 + t162 * t177, t162 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t179 ^ 2 / 0.2e1, t170 ^ 2 / 0.2e1, -t170 * t169, t170 * t176, -t169 * t176, t176 ^ 2 / 0.2e1, t159 * t169 + t198 * t176, t159 * t170 - t206 * t176, t164 ^ 2 / 0.2e1, -t164 * t163, t164 * t168, -t163 * t168, t168 ^ 2 / 0.2e1, (-t191 * t156 + t194 * t157) * t168 + t155 * t163, -(t194 * t156 + t191 * t157) * t168 + t155 * t164;];
T_reg = t1;
