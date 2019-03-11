% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:45
% EndTime: 2019-03-09 03:22:45
% DurationCPUTime: 0.12s
% Computational Cost: add. (228->40), mult. (465->83), div. (0->0), fcn. (274->6), ass. (0->33)
t145 = sin(pkin(9));
t146 = cos(pkin(9));
t148 = sin(qJ(3));
t149 = cos(qJ(3));
t137 = (t145 * t149 + t146 * t148) * qJD(1);
t150 = qJD(1) ^ 2;
t159 = t150 / 0.2e1;
t158 = cos(qJ(5));
t157 = t150 * qJ(2);
t140 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t152 = -qJ(4) * qJD(1) + t140;
t134 = qJD(3) * pkin(3) + t152 * t149;
t135 = t152 * t148;
t127 = t145 * t134 + t146 * t135;
t123 = qJD(3) * pkin(8) + t127;
t138 = (-t145 * t148 + t146 * t149) * qJD(1);
t139 = qJD(4) + (pkin(3) * t148 + qJ(2)) * qJD(1);
t128 = t137 * pkin(4) - t138 * pkin(8) + t139;
t147 = sin(qJ(5));
t156 = t158 * t123 + t147 * t128;
t155 = qJD(3) * t140;
t154 = qJD(1) * qJD(3);
t153 = -t147 * t123 + t158 * t128;
t126 = t146 * t134 - t145 * t135;
t122 = -qJD(3) * pkin(4) - t126;
t142 = -qJD(1) * pkin(1) + qJD(2);
t136 = qJD(5) + t137;
t130 = t147 * qJD(3) + t158 * t138;
t129 = -t158 * qJD(3) + t147 * t138;
t120 = t129 * pkin(5) + qJD(6) + t122;
t119 = -t129 * qJ(6) + t156;
t118 = t136 * pkin(5) - t130 * qJ(6) + t153;
t1 = [t159, 0, 0, t142 * qJD(1), t157, qJ(2) ^ 2 * t159 + t142 ^ 2 / 0.2e1, t149 ^ 2 * t159, -t149 * t150 * t148, t149 * t154, -t148 * t154, qJD(3) ^ 2 / 0.2e1, t148 * t157 + t149 * t155, -t148 * t155 + t149 * t157, -t126 * t138 - t127 * t137, t127 ^ 2 / 0.2e1 + t126 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1, t130 ^ 2 / 0.2e1, -t130 * t129, t130 * t136, -t129 * t136, t136 ^ 2 / 0.2e1, t122 * t129 + t153 * t136, t122 * t130 - t156 * t136, -t118 * t130 - t119 * t129, t119 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1;];
T_reg  = t1;
