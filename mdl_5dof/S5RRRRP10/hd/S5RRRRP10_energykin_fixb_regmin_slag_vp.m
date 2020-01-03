% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:20
% EndTime: 2019-12-31 22:12:20
% DurationCPUTime: 0.13s
% Computational Cost: add. (270->39), mult. (666->85), div. (0->0), fcn. (494->8), ass. (0->38)
t186 = cos(qJ(4));
t165 = sin(pkin(5));
t172 = qJD(1) ^ 2;
t185 = t165 ^ 2 * t172;
t180 = cos(pkin(5)) * qJD(1);
t163 = qJD(2) + t180;
t171 = cos(qJ(2));
t169 = sin(qJ(2));
t181 = qJD(1) * t165;
t177 = t169 * t181;
t179 = pkin(1) * t180;
t173 = -pkin(7) * t177 + t171 * t179;
t150 = -t163 * pkin(2) - t173;
t168 = sin(qJ(3));
t170 = cos(qJ(3));
t154 = -t170 * t163 + t168 * t177;
t155 = t168 * t163 + t170 * t177;
t141 = t154 * pkin(3) - t155 * pkin(9) + t150;
t176 = t171 * t181;
t158 = -qJD(3) + t176;
t182 = pkin(7) * t176 + t169 * t179;
t151 = t163 * pkin(8) + t182;
t152 = (-pkin(2) * t171 - pkin(8) * t169 - pkin(1)) * t181;
t183 = t170 * t151 + t168 * t152;
t144 = -t158 * pkin(9) + t183;
t167 = sin(qJ(4));
t184 = t167 * t141 + t186 * t144;
t178 = t171 * t185;
t175 = t186 * t141 - t167 * t144;
t174 = -t168 * t151 + t170 * t152;
t143 = t158 * pkin(3) - t174;
t153 = qJD(4) + t154;
t146 = t186 * t155 - t167 * t158;
t145 = t167 * t155 + t186 * t158;
t138 = t145 * pkin(4) + qJD(5) + t143;
t137 = -t145 * qJ(5) + t184;
t136 = t153 * pkin(4) - t146 * qJ(5) + t175;
t1 = [t172 / 0.2e1, 0, 0, t169 ^ 2 * t185 / 0.2e1, t169 * t178, t163 * t177, t163 * t176, t163 ^ 2 / 0.2e1, pkin(1) * t178 + t173 * t163, -pkin(1) * t169 * t185 - t182 * t163, t155 ^ 2 / 0.2e1, -t155 * t154, -t155 * t158, t154 * t158, t158 ^ 2 / 0.2e1, t150 * t154 - t174 * t158, t150 * t155 + t183 * t158, t146 ^ 2 / 0.2e1, -t146 * t145, t146 * t153, -t145 * t153, t153 ^ 2 / 0.2e1, t143 * t145 + t175 * t153, t143 * t146 - t184 * t153, -t136 * t146 - t137 * t145, t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1;];
T_reg = t1;
