% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:15
% EndTime: 2019-03-09 05:10:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (607->56), mult. (1516->112), div. (0->0), fcn. (1168->10), ass. (0->46)
t207 = cos(qJ(3));
t206 = pkin(7) + qJ(2);
t205 = cos(pkin(11));
t188 = qJD(3) + qJD(4);
t190 = sin(pkin(10));
t191 = cos(pkin(10));
t194 = sin(qJ(3));
t179 = (t190 * t207 + t191 * t194) * qJD(1);
t202 = qJD(1) * t190;
t180 = t206 * t202;
t201 = qJD(1) * t191;
t181 = t206 * t201;
t199 = -t180 * t207 - t194 * t181;
t166 = qJD(3) * pkin(3) - t179 * pkin(8) + t199;
t178 = t194 * t202 - t201 * t207;
t203 = -t194 * t180 + t207 * t181;
t169 = -t178 * pkin(8) + t203;
t193 = sin(qJ(4));
t196 = cos(qJ(4));
t204 = t193 * t166 + t196 * t169;
t156 = t188 * qJ(5) + t204;
t171 = t196 * t178 + t193 * t179;
t172 = -t193 * t178 + t196 * t179;
t182 = qJD(2) + (-pkin(2) * t191 - pkin(1)) * qJD(1);
t173 = t178 * pkin(3) + t182;
t159 = t171 * pkin(4) - t172 * qJ(5) + t173;
t189 = sin(pkin(11));
t152 = t205 * t156 + t189 * t159;
t151 = -t189 * t156 + t205 * t159;
t200 = t196 * t166 - t193 * t169;
t155 = -t188 * pkin(4) + qJD(5) - t200;
t197 = qJD(1) ^ 2;
t195 = cos(qJ(6));
t192 = sin(qJ(6));
t187 = t191 ^ 2;
t186 = t190 ^ 2;
t185 = -qJD(1) * pkin(1) + qJD(2);
t170 = qJD(6) + t171;
t168 = t172 * t205 + t189 * t188;
t167 = t189 * t172 - t188 * t205;
t161 = -t192 * t167 + t195 * t168;
t160 = t195 * t167 + t192 * t168;
t153 = t167 * pkin(5) + t155;
t150 = -t167 * pkin(9) + t152;
t149 = t171 * pkin(5) - t168 * pkin(9) + t151;
t1 = [t197 / 0.2e1, 0, 0, -t185 * t201, t185 * t202 (t186 + t187) * t197 * qJ(2), t185 ^ 2 / 0.2e1 + (t187 / 0.2e1 + t186 / 0.2e1) * qJ(2) ^ 2 * t197, t179 ^ 2 / 0.2e1, -t179 * t178, t179 * qJD(3), -t178 * qJD(3), qJD(3) ^ 2 / 0.2e1, qJD(3) * t199 + t182 * t178, -qJD(3) * t203 + t182 * t179, t172 ^ 2 / 0.2e1, -t172 * t171, t172 * t188, -t171 * t188, t188 ^ 2 / 0.2e1, t173 * t171 + t188 * t200, t173 * t172 - t188 * t204, t151 * t171 + t155 * t167, -t152 * t171 + t155 * t168, -t151 * t168 - t152 * t167, t152 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t170, -t160 * t170, t170 ^ 2 / 0.2e1 (t195 * t149 - t192 * t150) * t170 + t153 * t160 -(t192 * t149 + t195 * t150) * t170 + t153 * t161;];
T_reg  = t1;
