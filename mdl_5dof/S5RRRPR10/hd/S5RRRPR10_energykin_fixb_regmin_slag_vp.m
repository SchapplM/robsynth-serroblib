% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:50
% EndTime: 2019-12-31 21:29:50
% DurationCPUTime: 0.11s
% Computational Cost: add. (365->43), mult. (907->93), div. (0->0), fcn. (695->10), ass. (0->42)
t207 = cos(qJ(3));
t187 = sin(pkin(5));
t195 = qJD(1) ^ 2;
t206 = t187 ^ 2 * t195;
t202 = cos(pkin(5)) * qJD(1);
t184 = qJD(2) + t202;
t191 = sin(qJ(3));
t192 = sin(qJ(2));
t203 = qJD(1) * t187;
t199 = t192 * t203;
t176 = t191 * t184 + t207 * t199;
t194 = cos(qJ(2));
t198 = t194 * t203;
t179 = -qJD(3) + t198;
t201 = pkin(1) * t202;
t204 = pkin(7) * t198 + t192 * t201;
t172 = t184 * pkin(8) + t204;
t174 = (-pkin(2) * t194 - pkin(8) * t192 - pkin(1)) * t203;
t197 = -t191 * t172 + t207 * t174;
t159 = -t179 * pkin(3) - t176 * qJ(4) + t197;
t175 = -t207 * t184 + t191 * t199;
t205 = t207 * t172 + t191 * t174;
t161 = -t175 * qJ(4) + t205;
t186 = sin(pkin(10));
t188 = cos(pkin(10));
t156 = t186 * t159 + t188 * t161;
t200 = t194 * t206;
t166 = -t188 * t175 - t186 * t176;
t155 = t188 * t159 - t186 * t161;
t196 = -pkin(7) * t199 + t194 * t201;
t171 = -t184 * pkin(2) - t196;
t165 = t175 * pkin(3) + qJD(4) + t171;
t193 = cos(qJ(5));
t190 = sin(qJ(5));
t167 = -t186 * t175 + t188 * t176;
t164 = qJD(5) - t166;
t163 = t193 * t167 - t190 * t179;
t162 = t190 * t167 + t193 * t179;
t157 = -t166 * pkin(4) - t167 * pkin(9) + t165;
t154 = -t179 * pkin(9) + t156;
t153 = t179 * pkin(4) - t155;
t1 = [t195 / 0.2e1, 0, 0, t192 ^ 2 * t206 / 0.2e1, t192 * t200, t184 * t199, t184 * t198, t184 ^ 2 / 0.2e1, pkin(1) * t200 + t196 * t184, -pkin(1) * t192 * t206 - t204 * t184, t176 ^ 2 / 0.2e1, -t176 * t175, -t176 * t179, t175 * t179, t179 ^ 2 / 0.2e1, t171 * t175 - t197 * t179, t171 * t176 + t205 * t179, -t155 * t167 + t156 * t166, t156 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1, t163 ^ 2 / 0.2e1, -t163 * t162, t163 * t164, -t162 * t164, t164 ^ 2 / 0.2e1, (-t190 * t154 + t193 * t157) * t164 + t153 * t162, -(t193 * t154 + t190 * t157) * t164 + t153 * t163;];
T_reg = t1;
