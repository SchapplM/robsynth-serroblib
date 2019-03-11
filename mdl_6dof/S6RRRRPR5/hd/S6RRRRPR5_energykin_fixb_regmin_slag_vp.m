% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:16:16
% EndTime: 2019-03-09 22:16:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (439->51), mult. (924->105), div. (0->0), fcn. (657->8), ass. (0->46)
t208 = -pkin(4) - pkin(5);
t207 = -pkin(8) - pkin(7);
t191 = qJD(1) ^ 2;
t206 = t191 / 0.2e1;
t205 = cos(qJ(3));
t190 = cos(qJ(2));
t204 = t190 * t191;
t186 = sin(qJ(3));
t199 = qJD(1) * t190;
t187 = sin(qJ(2));
t200 = qJD(1) * t187;
t174 = t186 * t200 - t205 * t199;
t175 = (t186 * t190 + t205 * t187) * qJD(1);
t180 = (-pkin(2) * t190 - pkin(1)) * qJD(1);
t160 = t174 * pkin(3) - t175 * pkin(9) + t180;
t183 = qJD(2) + qJD(3);
t178 = qJD(2) * pkin(2) + t207 * t200;
t179 = t207 * t199;
t201 = t186 * t178 - t205 * t179;
t164 = t183 * pkin(9) + t201;
t185 = sin(qJ(4));
t189 = cos(qJ(4));
t203 = t185 * t160 + t189 * t164;
t202 = t205 * t178 + t186 * t179;
t198 = qJD(1) * qJD(2);
t173 = qJD(4) + t174;
t155 = t173 * qJ(5) + t203;
t197 = t187 * t198;
t196 = t190 * t198;
t195 = t189 * t160 - t185 * t164;
t194 = t183 * pkin(3) + t202;
t193 = qJD(5) - t195;
t167 = t189 * t175 + t185 * t183;
t192 = t167 * qJ(5) + t194;
t188 = cos(qJ(6));
t184 = sin(qJ(6));
t169 = -qJD(6) + t173;
t166 = t185 * t175 - t189 * t183;
t158 = t184 * t166 + t188 * t167;
t157 = -t188 * t166 + t184 * t167;
t156 = t166 * pkin(4) - t192;
t154 = -t173 * pkin(4) + t193;
t153 = t208 * t166 + t192;
t152 = t166 * pkin(10) + t155;
t151 = -t167 * pkin(10) + t208 * t173 + t193;
t1 = [t206, 0, 0, t187 ^ 2 * t206, t187 * t204, t197, t196, qJD(2) ^ 2 / 0.2e1, pkin(1) * t204 - pkin(7) * t197, -t191 * pkin(1) * t187 - pkin(7) * t196, t175 ^ 2 / 0.2e1, -t175 * t174, t175 * t183, -t174 * t183, t183 ^ 2 / 0.2e1, t180 * t174 + t202 * t183, t180 * t175 - t201 * t183, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t173, -t166 * t173, t173 ^ 2 / 0.2e1, -t166 * t194 + t195 * t173, -t167 * t194 - t203 * t173, -t154 * t173 + t156 * t166, t154 * t167 - t155 * t166, t155 * t173 - t156 * t167, t155 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t154 ^ 2 / 0.2e1, t158 ^ 2 / 0.2e1, -t157 * t158, -t169 * t158, t157 * t169, t169 ^ 2 / 0.2e1 -(t188 * t151 - t184 * t152) * t169 + t153 * t157 (t184 * t151 + t188 * t152) * t169 + t153 * t158;];
T_reg  = t1;
