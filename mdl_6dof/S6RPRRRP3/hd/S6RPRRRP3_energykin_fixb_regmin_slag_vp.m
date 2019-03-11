% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:15
% EndTime: 2019-03-09 06:05:15
% DurationCPUTime: 0.12s
% Computational Cost: add. (336->44), mult. (700->97), div. (0->0), fcn. (444->8), ass. (0->38)
t171 = qJD(1) ^ 2;
t184 = t171 / 0.2e1;
t183 = cos(qJ(4));
t167 = sin(qJ(4));
t168 = sin(qJ(3));
t179 = qJD(1) * t168;
t155 = t167 * qJD(3) + t183 * t179;
t170 = cos(qJ(3));
t178 = t170 * qJD(1);
t160 = -qJD(4) + t178;
t164 = sin(pkin(10));
t156 = (pkin(1) * t164 + pkin(7)) * qJD(1);
t180 = t168 * qJD(2) + t170 * t156;
t149 = qJD(3) * pkin(8) + t180;
t165 = cos(pkin(10));
t176 = -pkin(1) * t165 - pkin(2);
t150 = (-pkin(3) * t170 - pkin(8) * t168 + t176) * qJD(1);
t175 = -t167 * t149 + t183 * t150;
t139 = -t160 * pkin(4) - t155 * pkin(9) + t175;
t154 = -t183 * qJD(3) + t167 * t179;
t181 = t183 * t149 + t167 * t150;
t141 = -t154 * pkin(9) + t181;
t166 = sin(qJ(5));
t169 = cos(qJ(5));
t182 = t166 * t139 + t169 * t141;
t177 = qJD(1) * qJD(3);
t174 = t170 * qJD(2) - t168 * t156;
t173 = t169 * t139 - t166 * t141;
t148 = -qJD(3) * pkin(3) - t174;
t142 = t154 * pkin(4) + t148;
t158 = -qJD(5) + t160;
t157 = t176 * qJD(1);
t144 = -t166 * t154 + t169 * t155;
t143 = t169 * t154 + t166 * t155;
t137 = t143 * pkin(5) - t144 * qJ(6) + t142;
t136 = -t158 * qJ(6) + t182;
t135 = t158 * pkin(5) + qJD(6) - t173;
t1 = [t184, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t164 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t171, t168 ^ 2 * t184, t170 * t171 * t168, t168 * t177, t170 * t177, qJD(3) ^ 2 / 0.2e1, t174 * qJD(3) - t157 * t178, -t180 * qJD(3) + t157 * t179, t155 ^ 2 / 0.2e1, -t155 * t154, -t155 * t160, t154 * t160, t160 ^ 2 / 0.2e1, t148 * t154 - t175 * t160, t148 * t155 + t181 * t160, t144 ^ 2 / 0.2e1, -t144 * t143, -t144 * t158, t143 * t158, t158 ^ 2 / 0.2e1, t142 * t143 - t158 * t173, t142 * t144 + t182 * t158, t135 * t158 + t137 * t143, t135 * t144 - t136 * t143, -t136 * t158 - t137 * t144, t136 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1;];
T_reg  = t1;
