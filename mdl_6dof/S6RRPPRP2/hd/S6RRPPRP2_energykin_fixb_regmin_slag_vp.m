% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:31:56
% EndTime: 2019-03-09 08:31:56
% DurationCPUTime: 0.12s
% Computational Cost: add. (269->45), mult. (646->91), div. (0->0), fcn. (413->6), ass. (0->40)
t183 = pkin(3) + pkin(8);
t169 = qJD(1) ^ 2;
t182 = t169 / 0.2e1;
t181 = cos(qJ(5));
t180 = pkin(7) + qJ(3);
t168 = cos(qJ(2));
t179 = t168 * t169;
t164 = sin(pkin(9));
t165 = cos(pkin(9));
t176 = qJD(1) * t168;
t167 = sin(qJ(2));
t177 = qJD(1) * t167;
t155 = t164 * t177 - t165 * t176;
t156 = (t164 * t168 + t165 * t167) * qJD(1);
t161 = qJD(3) + (-pkin(2) * t168 - pkin(1)) * qJD(1);
t170 = -t156 * qJ(4) + t161;
t140 = t183 * t155 + t170;
t159 = qJD(2) * pkin(2) - t180 * t177;
t160 = t180 * t176;
t148 = t165 * t159 - t164 * t160;
t171 = qJD(4) - t148;
t143 = t156 * pkin(4) - t183 * qJD(2) + t171;
t166 = sin(qJ(5));
t178 = t181 * t140 + t166 * t143;
t149 = t164 * t159 + t165 * t160;
t175 = qJD(1) * qJD(2);
t147 = -qJD(2) * qJ(4) - t149;
t174 = t167 * t175;
t173 = t168 * t175;
t172 = -t166 * t140 + t181 * t143;
t144 = -t155 * pkin(4) - t147;
t154 = qJD(5) + t156;
t151 = t181 * qJD(2) + t166 * t155;
t150 = t166 * qJD(2) - t181 * t155;
t146 = -qJD(2) * pkin(3) + t171;
t145 = t155 * pkin(3) + t170;
t138 = t150 * pkin(5) + qJD(6) + t144;
t137 = -t150 * qJ(6) + t178;
t136 = t154 * pkin(5) - t151 * qJ(6) + t172;
t1 = [t182, 0, 0, t167 ^ 2 * t182, t167 * t179, t174, t173, qJD(2) ^ 2 / 0.2e1, pkin(1) * t179 - pkin(7) * t174, -t169 * pkin(1) * t167 - pkin(7) * t173, -t148 * t156 - t149 * t155, t149 ^ 2 / 0.2e1 + t148 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1, t146 * t156 + t147 * t155, t146 * qJD(2) - t145 * t155, -t147 * qJD(2) - t145 * t156, t145 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1, t151 ^ 2 / 0.2e1, -t151 * t150, t151 * t154, -t150 * t154, t154 ^ 2 / 0.2e1, t144 * t150 + t154 * t172, t144 * t151 - t178 * t154, -t136 * t151 - t137 * t150, t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1;];
T_reg  = t1;
