% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:15:03
% EndTime: 2019-03-09 21:15:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (473->49), mult. (977->99), div. (0->0), fcn. (657->6), ass. (0->40)
t175 = qJD(1) ^ 2;
t189 = t175 / 0.2e1;
t188 = pkin(4) + qJ(6);
t174 = cos(qJ(2));
t187 = t174 * t175;
t170 = sin(qJ(3));
t173 = cos(qJ(3));
t171 = sin(qJ(2));
t184 = qJD(1) * t171;
t159 = t170 * qJD(2) + t173 * t184;
t183 = t174 * qJD(1);
t165 = -qJD(3) + t183;
t157 = (-pkin(2) * t174 - pkin(8) * t171 - pkin(1)) * qJD(1);
t162 = pkin(7) * t183 + qJD(2) * pkin(8);
t178 = t173 * t157 - t170 * t162;
t147 = -t165 * pkin(3) - t159 * pkin(9) + t178;
t158 = -t173 * qJD(2) + t170 * t184;
t185 = t170 * t157 + t173 * t162;
t150 = -t158 * pkin(9) + t185;
t169 = sin(qJ(4));
t172 = cos(qJ(4));
t186 = t169 * t147 + t172 * t150;
t182 = qJD(1) * qJD(2);
t181 = t174 * t182;
t180 = t171 * t182;
t161 = -qJD(2) * pkin(2) + pkin(7) * t184;
t179 = t172 * t147 - t169 * t150;
t163 = -qJD(4) + t165;
t144 = t163 * qJ(5) - t186;
t177 = qJD(5) - t179;
t153 = t158 * pkin(3) + t161;
t152 = -t169 * t158 + t172 * t159;
t176 = -t152 * qJ(5) + t153;
t151 = t172 * t158 + t169 * t159;
t145 = t151 * pkin(4) + t176;
t143 = t163 * pkin(4) + t177;
t142 = t188 * t151 + t176;
t141 = -t151 * pkin(5) + qJD(6) - t144;
t140 = t152 * pkin(5) + t188 * t163 + t177;
t1 = [t189, 0, 0, t171 ^ 2 * t189, t171 * t187, t180, t181, qJD(2) ^ 2 / 0.2e1, pkin(1) * t187 - pkin(7) * t180, -t175 * pkin(1) * t171 - pkin(7) * t181, t159 ^ 2 / 0.2e1, -t159 * t158, -t159 * t165, t158 * t165, t165 ^ 2 / 0.2e1, t161 * t158 - t178 * t165, t161 * t159 + t185 * t165, t152 ^ 2 / 0.2e1, -t152 * t151, -t152 * t163, t151 * t163, t163 ^ 2 / 0.2e1, t153 * t151 - t179 * t163, t153 * t152 + t186 * t163, t143 * t152 + t144 * t151, -t143 * t163 - t145 * t151, t144 * t163 - t145 * t152, t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1, t140 * t152 - t141 * t151, -t141 * t163 - t142 * t152, t140 * t163 + t142 * t151, t142 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1;];
T_reg  = t1;
