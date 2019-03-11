% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPP5
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:09:07
% EndTime: 2019-03-09 21:09:07
% DurationCPUTime: 0.15s
% Computational Cost: add. (473->49), mult. (977->99), div. (0->0), fcn. (657->6), ass. (0->40)
t191 = pkin(4) + pkin(5);
t175 = qJD(1) ^ 2;
t190 = t175 / 0.2e1;
t189 = cos(qJ(4));
t174 = cos(qJ(2));
t188 = t174 * t175;
t171 = sin(qJ(3));
t173 = cos(qJ(3));
t172 = sin(qJ(2));
t185 = qJD(1) * t172;
t159 = t171 * qJD(2) + t173 * t185;
t184 = t174 * qJD(1);
t166 = -qJD(3) + t184;
t157 = (-pkin(2) * t174 - pkin(8) * t172 - pkin(1)) * qJD(1);
t163 = pkin(7) * t184 + qJD(2) * pkin(8);
t180 = t173 * t157 - t171 * t163;
t147 = -t166 * pkin(3) - t159 * pkin(9) + t180;
t158 = -t173 * qJD(2) + t171 * t185;
t186 = t171 * t157 + t173 * t163;
t150 = -t158 * pkin(9) + t186;
t170 = sin(qJ(4));
t187 = t170 * t147 + t189 * t150;
t183 = qJD(1) * qJD(2);
t164 = -qJD(4) + t166;
t144 = -t164 * qJ(5) + t187;
t182 = t172 * t183;
t181 = t174 * t183;
t162 = -qJD(2) * pkin(2) + pkin(7) * t185;
t179 = t189 * t147 - t170 * t150;
t178 = -t158 * pkin(3) - t162;
t177 = qJD(5) - t179;
t152 = -t170 * t158 + t189 * t159;
t176 = t152 * qJ(5) + t178;
t151 = t189 * t158 + t170 * t159;
t145 = t151 * pkin(4) - t176;
t143 = t164 * pkin(4) + t177;
t142 = -t191 * t151 + qJD(6) + t176;
t141 = t151 * qJ(6) + t144;
t140 = -t152 * qJ(6) + t191 * t164 + t177;
t1 = [t190, 0, 0, t172 ^ 2 * t190, t172 * t188, t182, t181, qJD(2) ^ 2 / 0.2e1, pkin(1) * t188 - pkin(7) * t182, -t175 * pkin(1) * t172 - pkin(7) * t181, t159 ^ 2 / 0.2e1, -t159 * t158, -t159 * t166, t158 * t166, t166 ^ 2 / 0.2e1, t162 * t158 - t180 * t166, t162 * t159 + t186 * t166, t152 ^ 2 / 0.2e1, -t152 * t151, -t152 * t164, t151 * t164, t164 ^ 2 / 0.2e1, -t151 * t178 - t164 * t179, -t152 * t178 + t187 * t164, t143 * t164 + t145 * t151, t143 * t152 - t144 * t151, -t144 * t164 - t145 * t152, t144 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1, t140 * t164 - t142 * t151, -t141 * t164 + t142 * t152, -t140 * t152 + t141 * t151, t141 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1;];
T_reg  = t1;
