% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:59:13
% EndTime: 2019-03-09 21:59:13
% DurationCPUTime: 0.17s
% Computational Cost: add. (690->54), mult. (1548->113), div. (0->0), fcn. (1171->10), ass. (0->49)
t211 = -pkin(8) - pkin(7);
t197 = qJD(1) ^ 2;
t210 = t197 / 0.2e1;
t209 = cos(qJ(3));
t208 = cos(pkin(11));
t196 = cos(qJ(2));
t207 = t196 * t197;
t188 = qJD(2) + qJD(3);
t187 = qJD(4) + t188;
t192 = sin(qJ(3));
t193 = sin(qJ(2));
t180 = (t192 * t196 + t193 * t209) * qJD(1);
t204 = qJD(1) * t193;
t182 = qJD(2) * pkin(2) + t204 * t211;
t203 = qJD(1) * t196;
t183 = t211 * t203;
t198 = t182 * t209 + t183 * t192;
t165 = pkin(3) * t188 - pkin(9) * t180 + t198;
t179 = t192 * t204 - t203 * t209;
t205 = t182 * t192 - t183 * t209;
t169 = -pkin(9) * t179 + t205;
t191 = sin(qJ(4));
t195 = cos(qJ(4));
t206 = t165 * t191 + t169 * t195;
t158 = qJ(5) * t187 + t206;
t173 = t179 * t195 + t180 * t191;
t174 = -t179 * t191 + t180 * t195;
t184 = (-pkin(2) * t196 - pkin(1)) * qJD(1);
t175 = t179 * pkin(3) + t184;
t163 = t173 * pkin(4) - t174 * qJ(5) + t175;
t189 = sin(pkin(11));
t154 = t158 * t208 + t163 * t189;
t202 = qJD(1) * qJD(2);
t201 = t193 * t202;
t200 = t196 * t202;
t153 = -t158 * t189 + t163 * t208;
t199 = t165 * t195 - t169 * t191;
t157 = -pkin(4) * t187 + qJD(5) - t199;
t194 = cos(qJ(6));
t190 = sin(qJ(6));
t172 = qJD(6) + t173;
t171 = t174 * t208 + t187 * t189;
t170 = t174 * t189 - t187 * t208;
t162 = -t170 * t190 + t171 * t194;
t161 = t170 * t194 + t171 * t190;
t155 = pkin(5) * t170 + t157;
t152 = -pkin(10) * t170 + t154;
t151 = pkin(5) * t173 - pkin(10) * t171 + t153;
t1 = [t210, 0, 0, t193 ^ 2 * t210, t193 * t207, t201, t200, qJD(2) ^ 2 / 0.2e1, pkin(1) * t207 - pkin(7) * t201, -pkin(1) * t193 * t197 - pkin(7) * t200, t180 ^ 2 / 0.2e1, -t180 * t179, t180 * t188, -t179 * t188, t188 ^ 2 / 0.2e1, t179 * t184 + t188 * t198, t180 * t184 - t188 * t205, t174 ^ 2 / 0.2e1, -t174 * t173, t174 * t187, -t173 * t187, t187 ^ 2 / 0.2e1, t173 * t175 + t187 * t199, t174 * t175 - t187 * t206, t153 * t173 + t157 * t170, -t154 * t173 + t157 * t171, -t153 * t171 - t154 * t170, t154 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * t172, -t161 * t172, t172 ^ 2 / 0.2e1 (t151 * t194 - t152 * t190) * t172 + t155 * t161 -(t151 * t190 + t152 * t194) * t172 + t155 * t162;];
T_reg  = t1;
