% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:19:57
% EndTime: 2019-03-09 04:19:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (184->47), mult. (356->95), div. (0->0), fcn. (186->6), ass. (0->34)
t155 = qJD(1) ^ 2;
t165 = t155 / 0.2e1;
t164 = cos(qJ(5));
t163 = t155 * qJ(2);
t143 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t154 = cos(qJ(3));
t131 = qJD(4) + (pkin(4) * qJD(1) - t143) * t154 + (-pkin(3) - pkin(8)) * qJD(3);
t152 = sin(qJ(3));
t160 = qJD(1) * t152;
t161 = pkin(3) * t160 + qJD(1) * qJ(2);
t133 = (pkin(8) * t152 - qJ(4) * t154) * qJD(1) + t161;
t151 = sin(qJ(5));
t162 = t151 * t131 + t164 * t133;
t138 = -qJD(3) * qJ(4) - t152 * t143;
t159 = qJD(1) * t154;
t158 = qJD(3) * t143;
t157 = qJD(1) * qJD(3);
t156 = t164 * t131 - t133 * t151;
t145 = qJD(5) + t159;
t134 = -pkin(4) * t160 - t138;
t153 = cos(qJ(6));
t150 = sin(qJ(6));
t147 = -qJD(1) * pkin(1) + qJD(2);
t142 = qJD(6) + t145;
t140 = qJD(3) * t164 + t151 * t160;
t139 = qJD(3) * t151 - t160 * t164;
t137 = -qJ(4) * t159 + t161;
t135 = -qJD(3) * pkin(3) - t154 * t143 + qJD(4);
t128 = -t139 * t150 + t140 * t153;
t127 = t153 * t139 + t140 * t150;
t126 = pkin(5) * t139 + t134;
t125 = -pkin(9) * t139 + t162;
t124 = pkin(5) * t145 - pkin(9) * t140 + t156;
t1 = [t165, 0, 0, t147 * qJD(1), t163, qJ(2) ^ 2 * t165 + t147 ^ 2 / 0.2e1, t154 ^ 2 * t165, -t154 * t155 * t152, t154 * t157, -t152 * t157, qJD(3) ^ 2 / 0.2e1, t152 * t163 + t154 * t158, -t152 * t158 + t154 * t163 (t135 * t154 + t138 * t152) * qJD(1), qJD(3) * t135 - t137 * t160, -qJD(3) * t138 - t137 * t159, t137 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1, t140 ^ 2 / 0.2e1, -t140 * t139, t140 * t145, -t139 * t145, t145 ^ 2 / 0.2e1, t134 * t139 + t145 * t156, t134 * t140 - t145 * t162, t128 ^ 2 / 0.2e1, -t128 * t127, t128 * t142, -t127 * t142, t142 ^ 2 / 0.2e1 (t124 * t153 - t125 * t150) * t142 + t126 * t127 -(t124 * t150 + t125 * t153) * t142 + t126 * t128;];
T_reg  = t1;
