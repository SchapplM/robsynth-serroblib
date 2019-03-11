% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:30
% EndTime: 2019-03-09 03:59:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (286->45), mult. (599->94), div. (0->0), fcn. (396->8), ass. (0->38)
t165 = sin(pkin(10));
t166 = cos(pkin(10));
t169 = sin(qJ(3));
t171 = cos(qJ(3));
t157 = (t165 * t171 + t166 * t169) * qJD(1);
t172 = qJD(1) ^ 2;
t181 = t172 / 0.2e1;
t180 = cos(qJ(5));
t179 = t172 * qJ(2);
t160 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t174 = -qJ(4) * qJD(1) + t160;
t153 = qJD(3) * pkin(3) + t174 * t171;
t155 = t174 * t169;
t145 = t165 * t153 + t166 * t155;
t141 = qJD(3) * pkin(8) + t145;
t158 = (-t165 * t169 + t166 * t171) * qJD(1);
t159 = qJD(4) + (pkin(3) * t169 + qJ(2)) * qJD(1);
t146 = t157 * pkin(4) - t158 * pkin(8) + t159;
t168 = sin(qJ(5));
t178 = t180 * t141 + t168 * t146;
t177 = qJD(3) * t160;
t176 = qJD(1) * qJD(3);
t175 = -t168 * t141 + t180 * t146;
t144 = t166 * t153 - t165 * t155;
t140 = -qJD(3) * pkin(4) - t144;
t156 = qJD(5) + t157;
t170 = cos(qJ(6));
t167 = sin(qJ(6));
t162 = -qJD(1) * pkin(1) + qJD(2);
t154 = qJD(6) + t156;
t149 = t168 * qJD(3) + t180 * t158;
t148 = -t180 * qJD(3) + t168 * t158;
t138 = -t167 * t148 + t170 * t149;
t137 = t170 * t148 + t167 * t149;
t136 = t148 * pkin(5) + t140;
t135 = -t148 * pkin(9) + t178;
t134 = t156 * pkin(5) - t149 * pkin(9) + t175;
t1 = [t181, 0, 0, t162 * qJD(1), t179, qJ(2) ^ 2 * t181 + t162 ^ 2 / 0.2e1, t171 ^ 2 * t181, -t171 * t172 * t169, t171 * t176, -t169 * t176, qJD(3) ^ 2 / 0.2e1, t169 * t179 + t171 * t177, -t169 * t177 + t171 * t179, -t144 * t158 - t145 * t157, t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1, t149 ^ 2 / 0.2e1, -t149 * t148, t149 * t156, -t148 * t156, t156 ^ 2 / 0.2e1, t140 * t148 + t175 * t156, t140 * t149 - t178 * t156, t138 ^ 2 / 0.2e1, -t138 * t137, t138 * t154, -t137 * t154, t154 ^ 2 / 0.2e1 (t170 * t134 - t167 * t135) * t154 + t136 * t137 -(t167 * t134 + t170 * t135) * t154 + t136 * t138;];
T_reg  = t1;
