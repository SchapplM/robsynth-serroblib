% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP2
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
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:03:09
% EndTime: 2019-03-10 01:03:10
% DurationCPUTime: 0.14s
% Computational Cost: add. (601->50), mult. (1354->105), div. (0->0), fcn. (1011->8), ass. (0->45)
t198 = -pkin(8) - pkin(7);
t183 = qJD(1) ^ 2;
t197 = t183 / 0.2e1;
t196 = cos(qJ(3));
t182 = cos(qJ(2));
t195 = t182 * t183;
t175 = qJD(2) + qJD(3);
t174 = qJD(4) + t175;
t178 = sin(qJ(3));
t179 = sin(qJ(2));
t167 = (t178 * t182 + t196 * t179) * qJD(1);
t191 = qJD(1) * t179;
t169 = qJD(2) * pkin(2) + t198 * t191;
t190 = qJD(1) * t182;
t170 = t198 * t190;
t185 = t196 * t169 + t178 * t170;
t153 = t175 * pkin(3) - t167 * pkin(9) + t185;
t166 = t178 * t191 - t196 * t190;
t192 = t178 * t169 - t196 * t170;
t156 = -t166 * pkin(9) + t192;
t177 = sin(qJ(4));
t181 = cos(qJ(4));
t193 = t177 * t153 + t181 * t156;
t149 = t174 * pkin(10) + t193;
t160 = t181 * t166 + t177 * t167;
t161 = -t177 * t166 + t181 * t167;
t171 = (-pkin(2) * t182 - pkin(1)) * qJD(1);
t162 = t166 * pkin(3) + t171;
t151 = t160 * pkin(4) - t161 * pkin(10) + t162;
t176 = sin(qJ(5));
t180 = cos(qJ(5));
t194 = t180 * t149 + t176 * t151;
t189 = qJD(1) * qJD(2);
t188 = t179 * t189;
t187 = t182 * t189;
t186 = t181 * t153 - t177 * t156;
t184 = -t176 * t149 + t180 * t151;
t148 = -t174 * pkin(4) - t186;
t159 = qJD(5) + t160;
t158 = t180 * t161 + t176 * t174;
t157 = t176 * t161 - t180 * t174;
t146 = t157 * pkin(5) - t158 * qJ(6) + t148;
t145 = t159 * qJ(6) + t194;
t144 = -t159 * pkin(5) + qJD(6) - t184;
t1 = [t197, 0, 0, t179 ^ 2 * t197, t179 * t195, t188, t187, qJD(2) ^ 2 / 0.2e1, pkin(1) * t195 - pkin(7) * t188, -t183 * pkin(1) * t179 - pkin(7) * t187, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t175, -t166 * t175, t175 ^ 2 / 0.2e1, t171 * t166 + t185 * t175, t171 * t167 - t192 * t175, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t174, -t160 * t174, t174 ^ 2 / 0.2e1, t162 * t160 + t186 * t174, t162 * t161 - t193 * t174, t158 ^ 2 / 0.2e1, -t158 * t157, t158 * t159, -t157 * t159, t159 ^ 2 / 0.2e1, t148 * t157 + t184 * t159, t148 * t158 - t194 * t159, -t144 * t159 + t146 * t157, t144 * t158 - t145 * t157, t145 * t159 - t146 * t158, t145 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1;];
T_reg  = t1;
