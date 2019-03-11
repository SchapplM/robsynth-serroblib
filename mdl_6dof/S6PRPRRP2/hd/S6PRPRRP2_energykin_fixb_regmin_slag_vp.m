% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:05
% EndTime: 2019-03-08 20:03:05
% DurationCPUTime: 0.11s
% Computational Cost: add. (204->37), mult. (455->83), div. (0->0), fcn. (315->10), ass. (0->36)
t185 = qJD(2) ^ 2;
t196 = t185 / 0.2e1;
t184 = cos(qJ(2));
t193 = qJD(1) * sin(pkin(6));
t168 = qJD(2) * pkin(2) + t184 * t193;
t176 = sin(pkin(11));
t178 = cos(pkin(11));
t181 = sin(qJ(2));
t189 = t181 * t193;
t164 = t176 * t168 + t178 * t189;
t162 = qJD(2) * pkin(8) + t164;
t172 = cos(pkin(6)) * qJD(1) + qJD(3);
t180 = sin(qJ(4));
t183 = cos(qJ(4));
t194 = t183 * t162 + t180 * t172;
t156 = qJD(4) * pkin(9) + t194;
t163 = t178 * t168 - t176 * t189;
t158 = (-pkin(4) * t183 - pkin(9) * t180 - pkin(3)) * qJD(2) - t163;
t179 = sin(qJ(5));
t182 = cos(qJ(5));
t195 = t182 * t156 + t179 * t158;
t192 = qJD(2) * t180;
t191 = t183 * qJD(2);
t190 = qJD(2) * qJD(4);
t188 = qJD(2) * t193;
t187 = -t180 * t162 + t183 * t172;
t186 = -t179 * t156 + t182 * t158;
t155 = -qJD(4) * pkin(4) - t187;
t173 = -qJD(5) + t191;
t167 = t179 * qJD(4) + t182 * t192;
t166 = -t182 * qJD(4) + t179 * t192;
t161 = -qJD(2) * pkin(3) - t163;
t153 = t166 * pkin(5) - t167 * qJ(6) + t155;
t152 = -t173 * qJ(6) + t195;
t151 = t173 * pkin(5) + qJD(6) - t186;
t1 = [qJD(1) ^ 2 / 0.2e1, t196, t184 * t188, -t181 * t188, t164 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1, t180 ^ 2 * t196, t180 * t185 * t183, t180 * t190, t183 * t190, qJD(4) ^ 2 / 0.2e1, t187 * qJD(4) - t161 * t191, -t194 * qJD(4) + t161 * t192, t167 ^ 2 / 0.2e1, -t167 * t166, -t167 * t173, t166 * t173, t173 ^ 2 / 0.2e1, t155 * t166 - t186 * t173, t155 * t167 + t195 * t173, t151 * t173 + t153 * t166, t151 * t167 - t152 * t166, -t152 * t173 - t153 * t167, t152 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1;];
T_reg  = t1;
