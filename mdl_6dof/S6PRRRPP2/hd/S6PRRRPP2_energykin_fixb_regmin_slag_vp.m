% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPP2
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:40
% EndTime: 2019-03-08 22:52:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (248->42), mult. (530->88), div. (0->0), fcn. (351->8), ass. (0->37)
t188 = pkin(4) + pkin(5);
t171 = qJD(2) ^ 2;
t187 = t171 / 0.2e1;
t186 = cos(qJ(4));
t168 = sin(qJ(2));
t182 = qJD(1) * sin(pkin(6));
t156 = qJD(2) * pkin(8) + t168 * t182;
t167 = sin(qJ(3));
t169 = cos(qJ(3));
t181 = qJD(1) * cos(pkin(6));
t183 = t169 * t156 + t167 * t181;
t149 = qJD(3) * pkin(9) + t183;
t170 = cos(qJ(2));
t177 = t170 * t182;
t151 = -t177 + (-pkin(3) * t169 - pkin(9) * t167 - pkin(2)) * qJD(2);
t166 = sin(qJ(4));
t185 = t186 * t149 + t166 * t151;
t184 = -t167 * t156 + t169 * t181;
t180 = qJD(2) * t167;
t179 = t169 * qJD(2);
t178 = qJD(2) * qJD(3);
t161 = -qJD(4) + t179;
t144 = -t161 * qJ(5) + t185;
t176 = qJD(2) * t182;
t175 = qJD(3) * pkin(3) + t184;
t174 = -t166 * t149 + t186 * t151;
t173 = qJD(5) - t174;
t155 = t166 * qJD(3) + t186 * t180;
t172 = t155 * qJ(5) + t175;
t157 = -qJD(2) * pkin(2) - t177;
t154 = -t186 * qJD(3) + t166 * t180;
t145 = t154 * pkin(4) - t172;
t143 = t161 * pkin(4) + t173;
t142 = -t188 * t154 + qJD(6) + t172;
t141 = t154 * qJ(6) + t144;
t140 = -t155 * qJ(6) + t188 * t161 + t173;
t1 = [qJD(1) ^ 2 / 0.2e1, t187, t170 * t176, -t168 * t176, t167 ^ 2 * t187, t167 * t171 * t169, t167 * t178, t169 * t178, qJD(3) ^ 2 / 0.2e1, t184 * qJD(3) - t157 * t179, -t183 * qJD(3) + t157 * t180, t155 ^ 2 / 0.2e1, -t155 * t154, -t155 * t161, t154 * t161, t161 ^ 2 / 0.2e1, -t154 * t175 - t174 * t161, -t155 * t175 + t185 * t161, t143 * t161 + t145 * t154, t143 * t155 - t144 * t154, -t144 * t161 - t145 * t155, t144 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1, t140 * t161 - t142 * t154, -t141 * t161 + t142 * t155, -t140 * t155 + t141 * t154, t141 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1;];
T_reg  = t1;
