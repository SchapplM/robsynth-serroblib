% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:30
% EndTime: 2019-12-31 19:24:30
% DurationCPUTime: 0.13s
% Computational Cost: add. (347->43), mult. (951->90), div. (0->0), fcn. (645->6), ass. (0->37)
t179 = sin(qJ(2));
t193 = qJD(1) * t179;
t196 = cos(pkin(5));
t170 = qJD(2) * pkin(2) + (-t196 * qJ(3) - pkin(7)) * t193;
t177 = sin(pkin(5));
t180 = cos(qJ(2));
t171 = (-qJ(3) * t177 * t179 - pkin(2) * t180 - pkin(1)) * qJD(1);
t200 = t196 * t170 + t171 * t177;
t187 = t180 * t196;
t191 = qJD(2) * t177;
t199 = qJD(1) * t187 + t191;
t192 = qJD(1) * t180;
t169 = pkin(7) * t192 + t199 * qJ(3);
t176 = sin(pkin(8));
t178 = cos(pkin(8));
t159 = -t176 * t169 + t200 * t178;
t181 = qJD(1) ^ 2;
t198 = t181 / 0.2e1;
t197 = pkin(3) + qJ(5);
t194 = t180 * t181;
t190 = qJD(1) * qJD(2);
t160 = t178 * t169 + t200 * t176;
t189 = t179 * t190;
t188 = t180 * t190;
t161 = -t177 * t170 + t196 * t171 + qJD(3);
t172 = -t196 * qJD(2) + t177 * t192;
t158 = t172 * qJ(4) - t160;
t163 = t176 * t191 + (t176 * t187 + t178 * t179) * qJD(1);
t183 = -t163 * qJ(4) + t161;
t182 = qJD(4) - t159;
t162 = t176 * t193 - t199 * t178;
t157 = t172 * pkin(3) + t182;
t156 = t162 * pkin(3) + t183;
t155 = -t162 * pkin(4) + qJD(5) - t158;
t154 = t197 * t162 + t183;
t153 = t163 * pkin(4) + t197 * t172 + t182;
t1 = [t198, 0, 0, t179 ^ 2 * t198, t179 * t194, t189, t188, qJD(2) ^ 2 / 0.2e1, pkin(1) * t194 - pkin(7) * t189, -t181 * pkin(1) * t179 - pkin(7) * t188, -t159 * t172 + t161 * t162, t160 * t172 + t161 * t163, -t159 * t163 - t160 * t162, t160 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1, t157 * t163 + t158 * t162, -t156 * t162 - t157 * t172, -t156 * t163 + t158 * t172, t156 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1, t153 * t163 - t155 * t162, -t154 * t163 - t155 * t172, t153 * t172 + t154 * t162, t154 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1;];
T_reg = t1;
