% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR3
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
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:03:59
% EndTime: 2019-03-09 22:03:59
% DurationCPUTime: 0.14s
% Computational Cost: add. (452->51), mult. (1040->105), div. (0->0), fcn. (758->8), ass. (0->46)
t192 = pkin(4) + pkin(10);
t191 = -pkin(8) - pkin(7);
t176 = qJD(1) ^ 2;
t190 = t176 / 0.2e1;
t189 = cos(qJ(3));
t175 = cos(qJ(2));
t188 = t175 * t176;
t171 = sin(qJ(3));
t172 = sin(qJ(2));
t161 = (t171 * t175 + t189 * t172) * qJD(1);
t168 = qJD(2) + qJD(3);
t185 = qJD(1) * t172;
t163 = qJD(2) * pkin(2) + t191 * t185;
t184 = qJD(1) * t175;
t164 = t191 * t184;
t179 = t189 * t163 + t171 * t164;
t146 = t168 * pkin(3) - t161 * pkin(9) + t179;
t160 = t171 * t185 - t189 * t184;
t186 = t171 * t163 - t189 * t164;
t149 = -t160 * pkin(9) + t186;
t170 = sin(qJ(4));
t174 = cos(qJ(4));
t187 = t170 * t146 + t174 * t149;
t183 = qJD(1) * qJD(2);
t182 = t172 * t183;
t181 = t175 * t183;
t180 = t174 * t146 - t170 * t149;
t167 = qJD(4) + t168;
t143 = -t167 * qJ(5) - t187;
t165 = (-pkin(2) * t175 - pkin(1)) * qJD(1);
t178 = qJD(5) - t180;
t155 = -t170 * t160 + t174 * t161;
t156 = t160 * pkin(3) + t165;
t177 = -t155 * qJ(5) + t156;
t173 = cos(qJ(6));
t169 = sin(qJ(6));
t154 = t174 * t160 + t170 * t161;
t153 = qJD(6) + t155;
t151 = t169 * t154 + t173 * t167;
t150 = -t173 * t154 + t169 * t167;
t144 = t154 * pkin(4) + t177;
t142 = -t167 * pkin(4) + t178;
t141 = t192 * t154 + t177;
t140 = -t154 * pkin(5) - t143;
t139 = t155 * pkin(5) - t192 * t167 + t178;
t1 = [t190, 0, 0, t172 ^ 2 * t190, t172 * t188, t182, t181, qJD(2) ^ 2 / 0.2e1, pkin(1) * t188 - pkin(7) * t182, -t176 * pkin(1) * t172 - pkin(7) * t181, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t168, -t160 * t168, t168 ^ 2 / 0.2e1, t165 * t160 + t179 * t168, t165 * t161 - t186 * t168, t155 ^ 2 / 0.2e1, -t155 * t154, t155 * t167, -t154 * t167, t167 ^ 2 / 0.2e1, t156 * t154 + t180 * t167, t156 * t155 - t187 * t167, t142 * t155 + t143 * t154, t142 * t167 - t144 * t154, -t143 * t167 - t144 * t155, t144 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1, t151 ^ 2 / 0.2e1, -t151 * t150, t151 * t153, -t150 * t153, t153 ^ 2 / 0.2e1 (t173 * t139 - t169 * t141) * t153 + t140 * t150 -(t169 * t139 + t173 * t141) * t153 + t140 * t151;];
T_reg  = t1;
