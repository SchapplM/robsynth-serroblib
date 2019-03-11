% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:54:54
% EndTime: 2019-03-09 21:54:54
% DurationCPUTime: 0.13s
% Computational Cost: add. (630->52), mult. (1519->109), div. (0->0), fcn. (1162->10), ass. (0->49)
t208 = -pkin(8) - pkin(7);
t195 = qJD(1) ^ 2;
t207 = t195 / 0.2e1;
t206 = cos(qJ(4));
t194 = cos(qJ(2));
t205 = t194 * t195;
t190 = sin(qJ(3));
t193 = cos(qJ(3));
t201 = qJD(1) * t194;
t191 = sin(qJ(2));
t202 = qJD(1) * t191;
t176 = t190 * t202 - t193 * t201;
t177 = (t190 * t194 + t191 * t193) * qJD(1);
t189 = sin(qJ(4));
t171 = -t189 * t176 + t206 * t177;
t185 = qJD(2) + qJD(3);
t184 = qJD(4) + t185;
t179 = qJD(2) * pkin(2) + t208 * t202;
t180 = t208 * t201;
t196 = t193 * t179 + t190 * t180;
t166 = t185 * pkin(3) - t177 * pkin(9) + t196;
t203 = t190 * t179 - t193 * t180;
t168 = -t176 * pkin(9) + t203;
t197 = t206 * t166 - t189 * t168;
t155 = t184 * pkin(4) - t171 * qJ(5) + t197;
t170 = t206 * t176 + t189 * t177;
t204 = t189 * t166 + t206 * t168;
t157 = -t170 * qJ(5) + t204;
t186 = sin(pkin(11));
t187 = cos(pkin(11));
t152 = t186 * t155 + t187 * t157;
t200 = qJD(1) * qJD(2);
t199 = t191 * t200;
t198 = t194 * t200;
t161 = -t187 * t170 - t186 * t171;
t181 = (-pkin(2) * t194 - pkin(1)) * qJD(1);
t151 = t187 * t155 - t186 * t157;
t172 = t176 * pkin(3) + t181;
t163 = t170 * pkin(4) + qJD(5) + t172;
t192 = cos(qJ(6));
t188 = sin(qJ(6));
t162 = -t186 * t170 + t187 * t171;
t160 = qJD(6) - t161;
t159 = t192 * t162 + t188 * t184;
t158 = t188 * t162 - t192 * t184;
t153 = -t161 * pkin(5) - t162 * pkin(10) + t163;
t150 = t184 * pkin(10) + t152;
t149 = -t184 * pkin(5) - t151;
t1 = [t207, 0, 0, t191 ^ 2 * t207, t191 * t205, t199, t198, qJD(2) ^ 2 / 0.2e1, pkin(1) * t205 - pkin(7) * t199, -t195 * pkin(1) * t191 - pkin(7) * t198, t177 ^ 2 / 0.2e1, -t177 * t176, t177 * t185, -t176 * t185, t185 ^ 2 / 0.2e1, t181 * t176 + t196 * t185, t181 * t177 - t203 * t185, t171 ^ 2 / 0.2e1, -t171 * t170, t171 * t184, -t170 * t184, t184 ^ 2 / 0.2e1, t172 * t170 + t197 * t184, t172 * t171 - t204 * t184, -t151 * t162 + t152 * t161, t152 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1, t159 ^ 2 / 0.2e1, -t159 * t158, t159 * t160, -t158 * t160, t160 ^ 2 / 0.2e1 (-t188 * t150 + t192 * t153) * t160 + t149 * t158 -(t192 * t150 + t188 * t153) * t160 + t149 * t159;];
T_reg  = t1;
