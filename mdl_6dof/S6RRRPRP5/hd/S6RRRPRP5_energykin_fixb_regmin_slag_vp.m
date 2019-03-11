% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:51:32
% EndTime: 2019-03-09 16:51:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (634->50), mult. (1371->103), div. (0->0), fcn. (978->8), ass. (0->43)
t190 = qJD(1) ^ 2;
t202 = t190 / 0.2e1;
t201 = cos(qJ(3));
t189 = cos(qJ(2));
t200 = t189 * t190;
t186 = sin(qJ(3));
t187 = sin(qJ(2));
t197 = qJD(1) * t187;
t173 = t186 * qJD(2) + t201 * t197;
t196 = t189 * qJD(1);
t179 = -qJD(3) + t196;
t171 = (-pkin(2) * t189 - pkin(8) * t187 - pkin(1)) * qJD(1);
t176 = pkin(7) * t196 + qJD(2) * pkin(8);
t192 = t201 * t171 - t186 * t176;
t162 = -t179 * pkin(3) - t173 * qJ(4) + t192;
t172 = -t201 * qJD(2) + t186 * t197;
t198 = t186 * t171 + t201 * t176;
t164 = -t172 * qJ(4) + t198;
t183 = sin(pkin(10));
t184 = cos(pkin(10));
t155 = t184 * t162 - t183 * t164;
t167 = -t183 * t172 + t184 * t173;
t152 = -t179 * pkin(4) - t167 * pkin(9) + t155;
t156 = t183 * t162 + t184 * t164;
t166 = -t184 * t172 - t183 * t173;
t154 = t166 * pkin(9) + t156;
t185 = sin(qJ(5));
t188 = cos(qJ(5));
t199 = t185 * t152 + t188 * t154;
t195 = qJD(1) * qJD(2);
t194 = t187 * t195;
t193 = t189 * t195;
t175 = -qJD(2) * pkin(2) + pkin(7) * t197;
t191 = t188 * t152 - t185 * t154;
t168 = t172 * pkin(3) + qJD(4) + t175;
t159 = -t166 * pkin(4) + t168;
t177 = -qJD(5) + t179;
t158 = t185 * t166 + t188 * t167;
t157 = -t188 * t166 + t185 * t167;
t150 = t157 * pkin(5) - t158 * qJ(6) + t159;
t149 = -t177 * qJ(6) + t199;
t148 = t177 * pkin(5) + qJD(6) - t191;
t1 = [t202, 0, 0, t187 ^ 2 * t202, t187 * t200, t194, t193, qJD(2) ^ 2 / 0.2e1, pkin(1) * t200 - pkin(7) * t194, -t190 * pkin(1) * t187 - pkin(7) * t193, t173 ^ 2 / 0.2e1, -t173 * t172, -t173 * t179, t172 * t179, t179 ^ 2 / 0.2e1, t175 * t172 - t192 * t179, t175 * t173 + t198 * t179, -t155 * t167 + t156 * t166, t156 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1 + t168 ^ 2 / 0.2e1, t158 ^ 2 / 0.2e1, -t158 * t157, -t158 * t177, t157 * t177, t177 ^ 2 / 0.2e1, t159 * t157 - t191 * t177, t159 * t158 + t199 * t177, t148 * t177 + t150 * t157, t148 * t158 - t149 * t157, -t149 * t177 - t150 * t158, t149 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1 + t148 ^ 2 / 0.2e1;];
T_reg  = t1;
