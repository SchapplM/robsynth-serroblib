% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:43:57
% EndTime: 2019-03-08 21:43:58
% DurationCPUTime: 0.11s
% Computational Cost: add. (156->41), mult. (354->86), div. (0->0), fcn. (216->8), ass. (0->37)
t184 = -pkin(3) - pkin(9);
t169 = qJD(2) ^ 2;
t183 = t169 / 0.2e1;
t165 = sin(qJ(2));
t180 = qJD(1) * sin(pkin(6));
t155 = qJD(2) * pkin(8) + t165 * t180;
t164 = sin(qJ(3));
t167 = cos(qJ(3));
t179 = qJD(1) * cos(pkin(6));
t171 = -t164 * t155 + t167 * t179;
t170 = qJD(4) - t171;
t177 = t164 * qJD(2);
t144 = pkin(4) * t177 + qJD(3) * t184 + t170;
t173 = -qJ(4) * t164 - pkin(2);
t168 = cos(qJ(2));
t175 = t168 * t180;
t147 = -t175 + (t167 * t184 + t173) * qJD(2);
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t182 = t163 * t144 + t166 * t147;
t181 = t167 * t155 + t164 * t179;
t178 = qJD(2) * t167;
t176 = qJD(2) * qJD(3);
t149 = -qJD(3) * qJ(4) - t181;
t174 = qJD(2) * t180;
t172 = t166 * t144 - t163 * t147;
t145 = pkin(4) * t178 - t149;
t158 = qJD(5) + t177;
t156 = -qJD(2) * pkin(2) - t175;
t154 = t166 * qJD(3) - t163 * t178;
t153 = t163 * qJD(3) + t166 * t178;
t150 = -t175 + (-pkin(3) * t167 + t173) * qJD(2);
t148 = -qJD(3) * pkin(3) + t170;
t141 = t153 * pkin(5) + qJD(6) + t145;
t140 = -t153 * qJ(6) + t182;
t139 = t158 * pkin(5) - t154 * qJ(6) + t172;
t1 = [qJD(1) ^ 2 / 0.2e1, t183, t168 * t174, -t165 * t174, t164 ^ 2 * t183, t164 * t169 * t167, t164 * t176, t167 * t176, qJD(3) ^ 2 / 0.2e1, qJD(3) * t171 - t156 * t178, -qJD(3) * t181 + t156 * t177 (t148 * t164 - t149 * t167) * qJD(2), t148 * qJD(3) + t150 * t178, -t149 * qJD(3) - t150 * t177, t150 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1 + t148 ^ 2 / 0.2e1, t154 ^ 2 / 0.2e1, -t154 * t153, t154 * t158, -t153 * t158, t158 ^ 2 / 0.2e1, t145 * t153 + t158 * t172, t145 * t154 - t158 * t182, -t139 * t154 - t140 * t153, t140 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1;];
T_reg  = t1;
