% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRP5
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
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:16
% EndTime: 2019-03-08 21:49:16
% DurationCPUTime: 0.11s
% Computational Cost: add. (204->43), mult. (438->90), div. (0->0), fcn. (269->8), ass. (0->37)
t185 = -pkin(3) - pkin(9);
t170 = qJD(2) ^ 2;
t184 = t170 / 0.2e1;
t166 = sin(qJ(2));
t181 = qJD(1) * sin(pkin(6));
t156 = qJD(2) * pkin(8) + t166 * t181;
t165 = sin(qJ(3));
t168 = cos(qJ(3));
t180 = qJD(1) * cos(pkin(6));
t173 = -t165 * t156 + t168 * t180;
t171 = qJD(4) - t173;
t178 = t165 * qJD(2);
t145 = pkin(4) * t178 + t185 * qJD(3) + t171;
t174 = -qJ(4) * t165 - pkin(2);
t169 = cos(qJ(2));
t176 = t169 * t181;
t148 = -t176 + (t185 * t168 + t174) * qJD(2);
t164 = sin(qJ(5));
t167 = cos(qJ(5));
t183 = t164 * t145 + t167 * t148;
t182 = t168 * t156 + t165 * t180;
t179 = qJD(2) * t168;
t177 = qJD(2) * qJD(3);
t150 = -qJD(3) * qJ(4) - t182;
t175 = qJD(2) * t181;
t146 = pkin(4) * t179 - t150;
t172 = t167 * t145 - t164 * t148;
t159 = qJD(5) + t178;
t157 = -qJD(2) * pkin(2) - t176;
t155 = t167 * qJD(3) - t164 * t179;
t154 = t164 * qJD(3) + t167 * t179;
t151 = -t176 + (-pkin(3) * t168 + t174) * qJD(2);
t149 = -qJD(3) * pkin(3) + t171;
t143 = t154 * pkin(5) - t155 * qJ(6) + t146;
t142 = t159 * qJ(6) + t183;
t141 = -t159 * pkin(5) + qJD(6) - t172;
t1 = [qJD(1) ^ 2 / 0.2e1, t184, t169 * t175, -t166 * t175, t165 ^ 2 * t184, t165 * t170 * t168, t165 * t177, t168 * t177, qJD(3) ^ 2 / 0.2e1, t173 * qJD(3) - t157 * t179, -t182 * qJD(3) + t157 * t178 (t149 * t165 - t150 * t168) * qJD(2), t149 * qJD(3) + t151 * t179, -t150 * qJD(3) - t151 * t178, t151 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1, t155 ^ 2 / 0.2e1, -t155 * t154, t155 * t159, -t154 * t159, t159 ^ 2 / 0.2e1, t146 * t154 + t159 * t172, t146 * t155 - t183 * t159, -t141 * t159 + t143 * t154, t141 * t155 - t142 * t154, t142 * t159 - t143 * t155, t142 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1;];
T_reg  = t1;
