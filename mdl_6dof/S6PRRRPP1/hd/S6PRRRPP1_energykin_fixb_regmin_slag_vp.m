% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:46:53
% EndTime: 2019-03-08 22:46:53
% DurationCPUTime: 0.14s
% Computational Cost: add. (344->43), mult. (750->92), div. (0->0), fcn. (529->10), ass. (0->40)
t199 = qJD(2) ^ 2;
t212 = t199 / 0.2e1;
t211 = cos(qJ(4));
t194 = sin(qJ(4));
t195 = sin(qJ(3));
t206 = qJD(2) * t195;
t183 = t194 * qJD(3) + t211 * t206;
t197 = cos(qJ(3));
t205 = t197 * qJD(2);
t187 = -qJD(4) + t205;
t196 = sin(qJ(2));
t208 = qJD(1) * sin(pkin(6));
t184 = qJD(2) * pkin(8) + t196 * t208;
t207 = qJD(1) * cos(pkin(6));
t209 = t197 * t184 + t195 * t207;
t175 = qJD(3) * pkin(9) + t209;
t198 = cos(qJ(2));
t203 = t198 * t208;
t178 = -t203 + (-pkin(3) * t197 - pkin(9) * t195 - pkin(2)) * qJD(2);
t201 = -t194 * t175 + t211 * t178;
t167 = -t187 * pkin(4) - t183 * qJ(5) + t201;
t182 = -t211 * qJD(3) + t194 * t206;
t210 = t211 * t175 + t194 * t178;
t169 = -t182 * qJ(5) + t210;
t190 = sin(pkin(11));
t192 = cos(pkin(11));
t164 = t190 * t167 + t192 * t169;
t204 = qJD(2) * qJD(3);
t202 = qJD(2) * t208;
t200 = -t195 * t184 + t197 * t207;
t163 = t192 * t167 - t190 * t169;
t174 = -qJD(3) * pkin(3) - t200;
t170 = t182 * pkin(4) + qJD(5) + t174;
t185 = -qJD(2) * pkin(2) - t203;
t172 = -t190 * t182 + t192 * t183;
t171 = t192 * t182 + t190 * t183;
t165 = t171 * pkin(5) - t172 * qJ(6) + t170;
t162 = -t187 * qJ(6) + t164;
t161 = pkin(5) * t187 + qJD(6) - t163;
t1 = [qJD(1) ^ 2 / 0.2e1, t212, t198 * t202, -t196 * t202, t195 ^ 2 * t212, t195 * t199 * t197, t195 * t204, t197 * t204, qJD(3) ^ 2 / 0.2e1, t200 * qJD(3) - t185 * t205, -t209 * qJD(3) + t185 * t206, t183 ^ 2 / 0.2e1, -t183 * t182, -t183 * t187, t182 * t187, t187 ^ 2 / 0.2e1, t174 * t182 - t201 * t187, t174 * t183 + t210 * t187, -t163 * t172 - t164 * t171, t164 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1 + t170 ^ 2 / 0.2e1, t161 * t187 + t165 * t171, t161 * t172 - t162 * t171, -t162 * t187 - t165 * t172, t162 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1;];
T_reg  = t1;
