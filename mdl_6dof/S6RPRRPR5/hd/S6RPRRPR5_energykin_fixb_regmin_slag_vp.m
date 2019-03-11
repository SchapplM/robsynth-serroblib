% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:33
% EndTime: 2019-03-09 05:13:33
% DurationCPUTime: 0.11s
% Computational Cost: add. (397->53), mult. (1018->104), div. (0->0), fcn. (755->8), ass. (0->43)
t188 = pkin(4) + pkin(9);
t187 = cos(qJ(3));
t186 = pkin(7) + qJ(2);
t169 = sin(pkin(10));
t170 = cos(pkin(10));
t173 = sin(qJ(3));
t160 = (t187 * t169 + t170 * t173) * qJD(1);
t183 = qJD(1) * t169;
t161 = t186 * t183;
t182 = qJD(1) * t170;
t162 = t186 * t182;
t180 = -t187 * t161 - t173 * t162;
t146 = qJD(3) * pkin(3) - t160 * pkin(8) + t180;
t159 = t173 * t183 - t187 * t182;
t184 = -t173 * t161 + t187 * t162;
t147 = -t159 * pkin(8) + t184;
t172 = sin(qJ(4));
t175 = cos(qJ(4));
t185 = t172 * t146 + t175 * t147;
t181 = t175 * t146 - t172 * t147;
t168 = qJD(3) + qJD(4);
t141 = -t168 * qJ(5) - t185;
t179 = qJD(5) - t181;
t153 = -t172 * t159 + t175 * t160;
t163 = qJD(2) + (-pkin(2) * t170 - pkin(1)) * qJD(1);
t154 = t159 * pkin(3) + t163;
t178 = -t153 * qJ(5) + t154;
t176 = qJD(1) ^ 2;
t174 = cos(qJ(6));
t171 = sin(qJ(6));
t167 = t170 ^ 2;
t166 = t169 ^ 2;
t165 = -qJD(1) * pkin(1) + qJD(2);
t152 = t175 * t159 + t172 * t160;
t151 = qJD(6) + t153;
t149 = t171 * t152 + t174 * t168;
t148 = -t174 * t152 + t171 * t168;
t142 = t152 * pkin(4) + t178;
t140 = -t168 * pkin(4) + t179;
t139 = t188 * t152 + t178;
t138 = -t152 * pkin(5) - t141;
t137 = t153 * pkin(5) - t188 * t168 + t179;
t1 = [t176 / 0.2e1, 0, 0, -t165 * t182, t165 * t183 (t166 + t167) * t176 * qJ(2), t165 ^ 2 / 0.2e1 + (t167 / 0.2e1 + t166 / 0.2e1) * qJ(2) ^ 2 * t176, t160 ^ 2 / 0.2e1, -t160 * t159, t160 * qJD(3), -t159 * qJD(3), qJD(3) ^ 2 / 0.2e1, t180 * qJD(3) + t163 * t159, -t184 * qJD(3) + t163 * t160, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t168, -t152 * t168, t168 ^ 2 / 0.2e1, t154 * t152 + t181 * t168, t154 * t153 - t185 * t168, t140 * t153 + t141 * t152, t140 * t168 - t142 * t152, -t141 * t168 - t142 * t153, t142 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1, t149 ^ 2 / 0.2e1, -t149 * t148, t149 * t151, -t148 * t151, t151 ^ 2 / 0.2e1 (t174 * t137 - t171 * t139) * t151 + t138 * t148 -(t171 * t137 + t174 * t139) * t151 + t138 * t149;];
T_reg  = t1;
