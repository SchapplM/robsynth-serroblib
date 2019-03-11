% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:29:50
% EndTime: 2019-03-09 04:29:50
% DurationCPUTime: 0.14s
% Computational Cost: add. (345->43), mult. (726->94), div. (0->0), fcn. (443->8), ass. (0->37)
t173 = qJD(1) ^ 2;
t184 = t173 / 0.2e1;
t183 = cos(qJ(4));
t170 = sin(qJ(4));
t171 = sin(qJ(3));
t180 = qJD(1) * t171;
t159 = t170 * qJD(3) + t183 * t180;
t172 = cos(qJ(3));
t179 = t172 * qJD(1);
t162 = -qJD(4) + t179;
t167 = sin(pkin(9));
t160 = (pkin(1) * t167 + pkin(7)) * qJD(1);
t181 = t171 * qJD(2) + t172 * t160;
t153 = qJD(3) * pkin(8) + t181;
t169 = cos(pkin(9));
t177 = -pkin(1) * t169 - pkin(2);
t154 = (-pkin(3) * t172 - pkin(8) * t171 + t177) * qJD(1);
t176 = -t170 * t153 + t183 * t154;
t143 = -t162 * pkin(4) - t159 * qJ(5) + t176;
t158 = -t183 * qJD(3) + t170 * t180;
t182 = t183 * t153 + t170 * t154;
t145 = -t158 * qJ(5) + t182;
t166 = sin(pkin(10));
t168 = cos(pkin(10));
t140 = t166 * t143 + t168 * t145;
t178 = qJD(1) * qJD(3);
t175 = t172 * qJD(2) - t171 * t160;
t139 = t168 * t143 - t166 * t145;
t152 = -qJD(3) * pkin(3) - t175;
t146 = t158 * pkin(4) + qJD(5) + t152;
t161 = t177 * qJD(1);
t148 = -t166 * t158 + t168 * t159;
t147 = t168 * t158 + t166 * t159;
t141 = t147 * pkin(5) - t148 * qJ(6) + t146;
t138 = -t162 * qJ(6) + t140;
t137 = t162 * pkin(5) + qJD(6) - t139;
t1 = [t184, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t167 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t173, t171 ^ 2 * t184, t172 * t173 * t171, t171 * t178, t172 * t178, qJD(3) ^ 2 / 0.2e1, t175 * qJD(3) - t161 * t179, -t181 * qJD(3) + t161 * t180, t159 ^ 2 / 0.2e1, -t159 * t158, -t159 * t162, t158 * t162, t162 ^ 2 / 0.2e1, t152 * t158 - t176 * t162, t152 * t159 + t182 * t162, -t139 * t148 - t140 * t147, t140 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1, t137 * t162 + t141 * t147, t137 * t148 - t138 * t147, -t138 * t162 - t141 * t148, t138 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1;];
T_reg  = t1;
