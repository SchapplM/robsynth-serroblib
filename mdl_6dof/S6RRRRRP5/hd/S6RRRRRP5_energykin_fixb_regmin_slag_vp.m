% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:23:03
% EndTime: 2019-03-10 01:23:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (500->49), mult. (1089->102), div. (0->0), fcn. (791->8), ass. (0->44)
t187 = qJD(1) ^ 2;
t202 = t187 / 0.2e1;
t201 = cos(qJ(3));
t200 = cos(qJ(5));
t186 = cos(qJ(2));
t199 = t186 * t187;
t183 = sin(qJ(3));
t184 = sin(qJ(2));
t195 = qJD(1) * t184;
t168 = -t201 * qJD(2) + t183 * t195;
t169 = t183 * qJD(2) + t201 * t195;
t182 = sin(qJ(4));
t185 = cos(qJ(4));
t162 = -t182 * t168 + t185 * t169;
t194 = t186 * qJD(1);
t177 = -qJD(3) + t194;
t175 = -qJD(4) + t177;
t167 = (-pkin(2) * t186 - pkin(8) * t184 - pkin(1)) * qJD(1);
t174 = pkin(7) * t194 + qJD(2) * pkin(8);
t188 = t201 * t167 - t183 * t174;
t157 = -t177 * pkin(3) - t169 * pkin(9) + t188;
t196 = t183 * t167 + t201 * t174;
t159 = -t168 * pkin(9) + t196;
t189 = t185 * t157 - t182 * t159;
t148 = -t175 * pkin(4) - t162 * pkin(10) + t189;
t161 = t185 * t168 + t182 * t169;
t197 = t182 * t157 + t185 * t159;
t150 = -t161 * pkin(10) + t197;
t181 = sin(qJ(5));
t198 = t181 * t148 + t200 * t150;
t193 = qJD(1) * qJD(2);
t192 = t184 * t193;
t191 = t186 * t193;
t173 = -qJD(2) * pkin(2) + pkin(7) * t195;
t190 = t200 * t148 - t181 * t150;
t163 = t168 * pkin(3) + t173;
t154 = t161 * pkin(4) + t163;
t171 = -qJD(5) + t175;
t153 = -t181 * t161 + t200 * t162;
t152 = t200 * t161 + t181 * t162;
t151 = t152 * pkin(5) + qJD(6) + t154;
t145 = -t152 * qJ(6) + t198;
t144 = -t171 * pkin(5) - t153 * qJ(6) + t190;
t1 = [t202, 0, 0, t184 ^ 2 * t202, t184 * t199, t192, t191, qJD(2) ^ 2 / 0.2e1, pkin(1) * t199 - pkin(7) * t192, -t187 * pkin(1) * t184 - pkin(7) * t191, t169 ^ 2 / 0.2e1, -t169 * t168, -t169 * t177, t168 * t177, t177 ^ 2 / 0.2e1, t173 * t168 - t188 * t177, t173 * t169 + t196 * t177, t162 ^ 2 / 0.2e1, -t162 * t161, -t162 * t175, t161 * t175, t175 ^ 2 / 0.2e1, t163 * t161 - t189 * t175, t163 * t162 + t197 * t175, t153 ^ 2 / 0.2e1, -t153 * t152, -t153 * t171, t152 * t171, t171 ^ 2 / 0.2e1, t154 * t152 - t190 * t171, t154 * t153 + t198 * t171, -t144 * t153 - t145 * t152, t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1;];
T_reg  = t1;
