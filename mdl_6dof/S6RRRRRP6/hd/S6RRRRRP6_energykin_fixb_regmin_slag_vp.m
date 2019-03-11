% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:31:26
% EndTime: 2019-03-10 01:31:26
% DurationCPUTime: 0.13s
% Computational Cost: add. (646->51), mult. (1362->106), div. (0->0), fcn. (989->8), ass. (0->44)
t192 = qJD(1) ^ 2;
t207 = t192 / 0.2e1;
t206 = cos(qJ(3));
t205 = cos(qJ(4));
t191 = cos(qJ(2));
t204 = t191 * t192;
t188 = sin(qJ(3));
t189 = sin(qJ(2));
t200 = qJD(1) * t189;
t173 = -t206 * qJD(2) + t188 * t200;
t174 = t188 * qJD(2) + t206 * t200;
t187 = sin(qJ(4));
t167 = -t187 * t173 + t205 * t174;
t199 = t191 * qJD(1);
t182 = -qJD(3) + t199;
t180 = -qJD(4) + t182;
t172 = (-pkin(2) * t191 - pkin(8) * t189 - pkin(1)) * qJD(1);
t179 = pkin(7) * t199 + qJD(2) * pkin(8);
t194 = t206 * t172 - t188 * t179;
t162 = -t182 * pkin(3) - t174 * pkin(9) + t194;
t201 = t188 * t172 + t206 * t179;
t164 = -t173 * pkin(9) + t201;
t195 = t205 * t162 - t187 * t164;
t154 = -t180 * pkin(4) - t167 * pkin(10) + t195;
t166 = t205 * t173 + t187 * t174;
t202 = t187 * t162 + t205 * t164;
t156 = -t166 * pkin(10) + t202;
t186 = sin(qJ(5));
t190 = cos(qJ(5));
t203 = t186 * t154 + t190 * t156;
t198 = qJD(1) * qJD(2);
t197 = t189 * t198;
t196 = t191 * t198;
t178 = -qJD(2) * pkin(2) + pkin(7) * t200;
t193 = t190 * t154 - t186 * t156;
t168 = t173 * pkin(3) + t178;
t159 = t166 * pkin(4) + t168;
t176 = -qJD(5) + t180;
t158 = -t186 * t166 + t190 * t167;
t157 = t190 * t166 + t186 * t167;
t152 = t157 * pkin(5) - t158 * qJ(6) + t159;
t151 = -t176 * qJ(6) + t203;
t150 = t176 * pkin(5) + qJD(6) - t193;
t1 = [t207, 0, 0, t189 ^ 2 * t207, t189 * t204, t197, t196, qJD(2) ^ 2 / 0.2e1, pkin(1) * t204 - pkin(7) * t197, -t192 * pkin(1) * t189 - pkin(7) * t196, t174 ^ 2 / 0.2e1, -t174 * t173, -t174 * t182, t173 * t182, t182 ^ 2 / 0.2e1, t178 * t173 - t194 * t182, t178 * t174 + t201 * t182, t167 ^ 2 / 0.2e1, -t167 * t166, -t167 * t180, t166 * t180, t180 ^ 2 / 0.2e1, t168 * t166 - t195 * t180, t168 * t167 + t202 * t180, t158 ^ 2 / 0.2e1, -t158 * t157, -t158 * t176, t157 * t176, t176 ^ 2 / 0.2e1, t159 * t157 - t193 * t176, t159 * t158 + t203 * t176, t150 * t176 + t152 * t157, t150 * t158 - t151 * t157, -t151 * t176 - t152 * t158, t151 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1;];
T_reg  = t1;
