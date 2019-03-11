% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:02:34
% EndTime: 2019-03-09 21:02:34
% DurationCPUTime: 0.12s
% Computational Cost: add. (670->50), mult. (1420->103), div. (0->0), fcn. (996->8), ass. (0->43)
t195 = qJD(1) ^ 2;
t208 = t195 / 0.2e1;
t207 = cos(qJ(3));
t206 = cos(qJ(4));
t194 = cos(qJ(2));
t205 = t194 * t195;
t192 = sin(qJ(3));
t193 = sin(qJ(2));
t202 = qJD(1) * t193;
t178 = -t207 * qJD(2) + t192 * t202;
t179 = t192 * qJD(2) + t207 * t202;
t191 = sin(qJ(4));
t172 = -t191 * t178 + t206 * t179;
t201 = t194 * qJD(1);
t185 = -qJD(3) + t201;
t183 = -qJD(4) + t185;
t177 = (-pkin(2) * t194 - pkin(8) * t193 - pkin(1)) * qJD(1);
t182 = pkin(7) * t201 + qJD(2) * pkin(8);
t196 = t207 * t177 - t192 * t182;
t167 = -t185 * pkin(3) - t179 * pkin(9) + t196;
t203 = t192 * t177 + t207 * t182;
t169 = -t178 * pkin(9) + t203;
t197 = t206 * t167 - t191 * t169;
t159 = -t183 * pkin(4) - t172 * qJ(5) + t197;
t171 = t206 * t178 + t191 * t179;
t204 = t191 * t167 + t206 * t169;
t161 = -t171 * qJ(5) + t204;
t189 = sin(pkin(10));
t190 = cos(pkin(10));
t156 = t189 * t159 + t190 * t161;
t200 = qJD(1) * qJD(2);
t199 = t193 * t200;
t198 = t194 * t200;
t181 = -qJD(2) * pkin(2) + pkin(7) * t202;
t155 = t190 * t159 - t189 * t161;
t173 = t178 * pkin(3) + t181;
t164 = t171 * pkin(4) + qJD(5) + t173;
t163 = -t189 * t171 + t190 * t172;
t162 = t190 * t171 + t189 * t172;
t157 = t162 * pkin(5) - t163 * qJ(6) + t164;
t154 = -t183 * qJ(6) + t156;
t153 = t183 * pkin(5) + qJD(6) - t155;
t1 = [t208, 0, 0, t193 ^ 2 * t208, t193 * t205, t199, t198, qJD(2) ^ 2 / 0.2e1, pkin(1) * t205 - pkin(7) * t199, -t195 * pkin(1) * t193 - pkin(7) * t198, t179 ^ 2 / 0.2e1, -t179 * t178, -t179 * t185, t178 * t185, t185 ^ 2 / 0.2e1, t181 * t178 - t196 * t185, t181 * t179 + t203 * t185, t172 ^ 2 / 0.2e1, -t172 * t171, -t172 * t183, t171 * t183, t183 ^ 2 / 0.2e1, t173 * t171 - t197 * t183, t173 * t172 + t204 * t183, -t155 * t163 - t156 * t162, t156 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1, t153 * t183 + t157 * t162, t153 * t163 - t154 * t162, -t154 * t183 - t157 * t163, t154 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1;];
T_reg  = t1;
