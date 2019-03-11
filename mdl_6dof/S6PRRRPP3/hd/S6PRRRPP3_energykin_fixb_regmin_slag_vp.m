% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:29
% EndTime: 2019-03-08 22:58:29
% DurationCPUTime: 0.11s
% Computational Cost: add. (248->42), mult. (530->88), div. (0->0), fcn. (351->8), ass. (0->37)
t170 = qJD(2) ^ 2;
t185 = t170 / 0.2e1;
t184 = pkin(4) + qJ(6);
t166 = sin(qJ(2));
t181 = qJD(1) * sin(pkin(6));
t156 = qJD(2) * pkin(8) + t166 * t181;
t165 = sin(qJ(3));
t168 = cos(qJ(3));
t180 = qJD(1) * cos(pkin(6));
t182 = t168 * t156 + t165 * t180;
t149 = qJD(3) * pkin(9) + t182;
t169 = cos(qJ(2));
t176 = t169 * t181;
t151 = -t176 + (-pkin(3) * t168 - pkin(9) * t165 - pkin(2)) * qJD(2);
t164 = sin(qJ(4));
t167 = cos(qJ(4));
t183 = t167 * t149 + t164 * t151;
t179 = qJD(2) * t165;
t178 = t168 * qJD(2);
t177 = qJD(2) * qJD(3);
t175 = qJD(2) * t181;
t174 = -t164 * t149 + t167 * t151;
t159 = -qJD(4) + t178;
t144 = t159 * qJ(5) - t183;
t173 = -t165 * t156 + t168 * t180;
t172 = qJD(5) - t174;
t148 = -qJD(3) * pkin(3) - t173;
t155 = t164 * qJD(3) + t167 * t179;
t171 = -t155 * qJ(5) + t148;
t157 = -qJD(2) * pkin(2) - t176;
t154 = -t167 * qJD(3) + t164 * t179;
t145 = t154 * pkin(4) + t171;
t143 = t159 * pkin(4) + t172;
t142 = t184 * t154 + t171;
t141 = -t154 * pkin(5) + qJD(6) - t144;
t140 = t155 * pkin(5) + t184 * t159 + t172;
t1 = [qJD(1) ^ 2 / 0.2e1, t185, t169 * t175, -t166 * t175, t165 ^ 2 * t185, t165 * t170 * t168, t165 * t177, t168 * t177, qJD(3) ^ 2 / 0.2e1, t173 * qJD(3) - t157 * t178, -t182 * qJD(3) + t157 * t179, t155 ^ 2 / 0.2e1, -t155 * t154, -t155 * t159, t154 * t159, t159 ^ 2 / 0.2e1, t148 * t154 - t174 * t159, t148 * t155 + t183 * t159, t143 * t155 + t144 * t154, -t143 * t159 - t145 * t154, t144 * t159 - t145 * t155, t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1, t140 * t155 - t141 * t154, -t141 * t159 - t142 * t155, t140 * t159 + t142 * t154, t142 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1;];
T_reg  = t1;
