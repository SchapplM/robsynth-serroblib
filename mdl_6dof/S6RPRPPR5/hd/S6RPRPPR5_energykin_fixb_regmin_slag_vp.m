% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:37
% EndTime: 2019-03-09 02:51:37
% DurationCPUTime: 0.13s
% Computational Cost: add. (371->54), mult. (908->105), div. (0->0), fcn. (625->8), ass. (0->42)
t195 = pkin(3) + qJ(5);
t194 = pkin(7) + qJ(2);
t193 = cos(pkin(10));
t182 = sin(qJ(3));
t184 = cos(qJ(3));
t180 = cos(pkin(9));
t190 = qJD(1) * t180;
t179 = sin(pkin(9));
t191 = qJD(1) * t179;
t168 = t182 * t191 - t184 * t190;
t169 = (t179 * t184 + t180 * t182) * qJD(1);
t172 = qJD(2) + (-pkin(2) * t180 - pkin(1)) * qJD(1);
t187 = -t169 * qJ(4) + t172;
t150 = t195 * t168 + t187;
t170 = t194 * t191;
t171 = t194 * t190;
t189 = -t184 * t170 - t182 * t171;
t188 = qJD(4) - t189;
t153 = t169 * pkin(4) - t195 * qJD(3) + t188;
t178 = sin(pkin(10));
t147 = t193 * t150 + t178 * t153;
t192 = -t182 * t170 + t184 * t171;
t159 = -qJD(3) * qJ(4) - t192;
t146 = -t178 * t150 + t193 * t153;
t156 = -t168 * pkin(4) + qJD(5) - t159;
t185 = qJD(1) ^ 2;
t183 = cos(qJ(6));
t181 = sin(qJ(6));
t176 = t180 ^ 2;
t175 = t179 ^ 2;
t174 = -qJD(1) * pkin(1) + qJD(2);
t164 = qJD(6) + t169;
t162 = t193 * qJD(3) + t178 * t168;
t161 = t178 * qJD(3) - t193 * t168;
t158 = -qJD(3) * pkin(3) + t188;
t157 = t168 * pkin(3) + t187;
t155 = -t181 * t161 + t183 * t162;
t154 = t183 * t161 + t181 * t162;
t148 = t161 * pkin(5) + t156;
t145 = -t161 * pkin(8) + t147;
t144 = t169 * pkin(5) - t162 * pkin(8) + t146;
t1 = [t185 / 0.2e1, 0, 0, -t174 * t190, t174 * t191 (t175 + t176) * t185 * qJ(2), t174 ^ 2 / 0.2e1 + (t176 / 0.2e1 + t175 / 0.2e1) * qJ(2) ^ 2 * t185, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * qJD(3), -t168 * qJD(3), qJD(3) ^ 2 / 0.2e1, t189 * qJD(3) + t172 * t168, -t192 * qJD(3) + t172 * t169, t158 * t169 + t159 * t168, t158 * qJD(3) - t157 * t168, -t159 * qJD(3) - t157 * t169, t157 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, t146 * t169 + t156 * t161, -t147 * t169 + t156 * t162, -t146 * t162 - t147 * t161, t147 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1, t155 ^ 2 / 0.2e1, -t155 * t154, t155 * t164, -t154 * t164, t164 ^ 2 / 0.2e1 (t183 * t144 - t181 * t145) * t164 + t148 * t154 -(t181 * t144 + t183 * t145) * t164 + t148 * t155;];
T_reg  = t1;
