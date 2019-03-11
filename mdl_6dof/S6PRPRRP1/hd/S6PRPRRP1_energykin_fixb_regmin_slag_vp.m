% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:28
% EndTime: 2019-03-08 19:58:28
% DurationCPUTime: 0.12s
% Computational Cost: add. (150->35), mult. (355->79), div. (0->0), fcn. (244->10), ass. (0->36)
t176 = qJD(2) ^ 2;
t188 = t176 / 0.2e1;
t187 = cos(qJ(5));
t175 = cos(qJ(2));
t184 = qJD(1) * sin(pkin(6));
t160 = qJD(2) * pkin(2) + t175 * t184;
t168 = sin(pkin(11));
t170 = cos(pkin(11));
t173 = sin(qJ(2));
t180 = t173 * t184;
t156 = t168 * t160 + t170 * t180;
t154 = qJD(2) * pkin(8) + t156;
t164 = cos(pkin(6)) * qJD(1) + qJD(3);
t172 = sin(qJ(4));
t174 = cos(qJ(4));
t185 = t174 * t154 + t172 * t164;
t147 = qJD(4) * pkin(9) + t185;
t155 = t170 * t160 - t168 * t180;
t150 = (-pkin(4) * t174 - pkin(9) * t172 - pkin(3)) * qJD(2) - t155;
t171 = sin(qJ(5));
t186 = t187 * t147 + t171 * t150;
t183 = qJD(2) * t172;
t182 = t174 * qJD(2);
t181 = qJD(2) * qJD(4);
t179 = qJD(2) * t184;
t178 = -t171 * t147 + t187 * t150;
t177 = -t172 * t154 + t174 * t164;
t146 = -qJD(4) * pkin(4) - t177;
t165 = -qJD(5) + t182;
t159 = t171 * qJD(4) + t183 * t187;
t158 = -qJD(4) * t187 + t171 * t183;
t153 = -qJD(2) * pkin(3) - t155;
t144 = t158 * pkin(5) + qJD(6) + t146;
t143 = -t158 * qJ(6) + t186;
t142 = -t165 * pkin(5) - t159 * qJ(6) + t178;
t1 = [qJD(1) ^ 2 / 0.2e1, t188, t175 * t179, -t173 * t179, t156 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1, t172 ^ 2 * t188, t172 * t176 * t174, t172 * t181, t174 * t181, qJD(4) ^ 2 / 0.2e1, qJD(4) * t177 - t153 * t182, -qJD(4) * t185 + t153 * t183, t159 ^ 2 / 0.2e1, -t159 * t158, -t159 * t165, t158 * t165, t165 ^ 2 / 0.2e1, t146 * t158 - t165 * t178, t146 * t159 + t165 * t186, -t142 * t159 - t143 * t158, t143 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1;];
T_reg  = t1;
