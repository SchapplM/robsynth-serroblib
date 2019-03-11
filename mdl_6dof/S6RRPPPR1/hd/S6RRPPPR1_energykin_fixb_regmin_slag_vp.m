% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:17
% EndTime: 2019-03-09 08:08:17
% DurationCPUTime: 0.13s
% Computational Cost: add. (421->51), mult. (1011->103), div. (0->0), fcn. (692->8), ass. (0->44)
t216 = -pkin(4) - pkin(5);
t203 = qJD(1) ^ 2;
t215 = t203 / 0.2e1;
t214 = pkin(7) + qJ(3);
t213 = cos(pkin(9));
t202 = cos(qJ(2));
t212 = t202 * t203;
t197 = sin(pkin(9));
t210 = qJD(1) * t202;
t200 = sin(qJ(2));
t211 = qJD(1) * t200;
t187 = t197 * t211 - t213 * t210;
t188 = (t197 * t202 + t213 * t200) * qJD(1);
t193 = qJD(3) + (-pkin(2) * t202 - pkin(1)) * qJD(1);
t172 = t187 * pkin(3) - t188 * qJ(4) + t193;
t191 = qJD(2) * pkin(2) - t214 * t211;
t192 = t214 * t210;
t178 = t197 * t191 + t213 * t192;
t176 = qJD(2) * qJ(4) + t178;
t196 = sin(pkin(10));
t198 = cos(pkin(10));
t168 = t196 * t172 + t198 * t176;
t177 = t213 * t191 - t197 * t192;
t209 = qJD(1) * qJD(2);
t165 = t187 * qJ(5) + t168;
t208 = t200 * t209;
t207 = t202 * t209;
t167 = t198 * t172 - t196 * t176;
t206 = qJD(5) - t167;
t205 = qJD(2) * pkin(3) - qJD(4) + t177;
t181 = t196 * qJD(2) + t198 * t188;
t204 = t181 * qJ(5) + t205;
t201 = cos(qJ(6));
t199 = sin(qJ(6));
t185 = -qJD(6) + t187;
t180 = -t198 * qJD(2) + t196 * t188;
t170 = t199 * t180 + t201 * t181;
t169 = -t201 * t180 + t199 * t181;
t166 = t180 * pkin(4) - t204;
t164 = -t187 * pkin(4) + t206;
t163 = t216 * t180 + t204;
t162 = t180 * pkin(8) + t165;
t161 = -t181 * pkin(8) + t216 * t187 + t206;
t1 = [t215, 0, 0, t200 ^ 2 * t215, t200 * t212, t208, t207, qJD(2) ^ 2 / 0.2e1, pkin(1) * t212 - pkin(7) * t208, -t203 * pkin(1) * t200 - pkin(7) * t207, -t177 * t188 - t178 * t187, t178 ^ 2 / 0.2e1 + t177 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1, t167 * t187 - t180 * t205, -t168 * t187 - t181 * t205, -t167 * t181 - t168 * t180, t168 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1 + t205 ^ 2 / 0.2e1, -t164 * t187 + t166 * t180, t164 * t181 - t165 * t180, t165 * t187 - t166 * t181, t165 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1, t170 ^ 2 / 0.2e1, -t170 * t169, -t170 * t185, t169 * t185, t185 ^ 2 / 0.2e1 -(t201 * t161 - t199 * t162) * t185 + t163 * t169 (t199 * t161 + t201 * t162) * t185 + t163 * t170;];
T_reg  = t1;
