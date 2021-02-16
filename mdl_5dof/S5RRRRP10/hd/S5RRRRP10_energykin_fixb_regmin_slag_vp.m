% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:35:13
% EndTime: 2021-01-16 00:35:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (364->41), mult. (872->89), div. (0->0), fcn. (653->8), ass. (0->38)
t194 = cos(qJ(4));
t173 = sin(pkin(5));
t180 = qJD(1) ^ 2;
t193 = t173 ^ 2 * t180;
t188 = cos(pkin(5)) * qJD(1);
t171 = qJD(2) + t188;
t179 = cos(qJ(2));
t177 = sin(qJ(2));
t189 = qJD(1) * t173;
t185 = t177 * t189;
t187 = pkin(1) * t188;
t181 = -pkin(7) * t185 + t179 * t187;
t158 = -t171 * pkin(2) - t181;
t176 = sin(qJ(3));
t178 = cos(qJ(3));
t162 = -t178 * t171 + t176 * t185;
t163 = t176 * t171 + t178 * t185;
t149 = t162 * pkin(3) - t163 * pkin(9) + t158;
t184 = t179 * t189;
t166 = -qJD(3) + t184;
t190 = pkin(7) * t184 + t177 * t187;
t159 = t171 * pkin(8) + t190;
t160 = (-pkin(2) * t179 - pkin(8) * t177 - pkin(1)) * t189;
t191 = t178 * t159 + t176 * t160;
t152 = -t166 * pkin(9) + t191;
t175 = sin(qJ(4));
t192 = t175 * t149 + t194 * t152;
t186 = t179 * t193;
t183 = t194 * t149 - t175 * t152;
t182 = -t176 * t159 + t178 * t160;
t151 = t166 * pkin(3) - t182;
t161 = qJD(4) + t162;
t154 = t194 * t163 - t175 * t166;
t153 = t175 * t163 + t194 * t166;
t146 = t153 * pkin(4) + qJD(5) + t151;
t145 = -t153 * qJ(5) + t192;
t144 = t161 * pkin(4) - t154 * qJ(5) + t183;
t1 = [t180 / 0.2e1, 0, 0, t177 ^ 2 * t193 / 0.2e1, t177 * t186, t171 * t185, t171 * t184, t171 ^ 2 / 0.2e1, pkin(1) * t186 + t181 * t171, -pkin(1) * t177 * t193 - t190 * t171, t163 ^ 2 / 0.2e1, -t163 * t162, -t163 * t166, t162 * t166, t166 ^ 2 / 0.2e1, t158 * t162 - t182 * t166, t158 * t163 + t191 * t166, t154 ^ 2 / 0.2e1, -t154 * t153, t154 * t161, -t153 * t161, t161 ^ 2 / 0.2e1, t151 * t153 + t183 * t161, t151 * t154 - t192 * t161, t144 * t161 + t146 * t153, -t145 * t161 + t146 * t154, -t144 * t154 - t145 * t153, t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1;];
T_reg = t1;
