% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:36:57
% EndTime: 2019-03-09 04:36:58
% DurationCPUTime: 0.08s
% Computational Cost: add. (249->42), mult. (512->90), div. (0->0), fcn. (285->6), ass. (0->34)
t148 = qJD(1) ^ 2;
t161 = t148 / 0.2e1;
t160 = pkin(4) + qJ(6);
t142 = sin(pkin(9));
t136 = (pkin(1) * t142 + pkin(7)) * qJD(1);
t145 = sin(qJ(3));
t147 = cos(qJ(3));
t158 = t145 * qJD(2) + t147 * t136;
t130 = qJD(3) * pkin(8) + t158;
t143 = cos(pkin(9));
t154 = -pkin(1) * t143 - pkin(2);
t131 = (-pkin(3) * t147 - pkin(8) * t145 + t154) * qJD(1);
t144 = sin(qJ(4));
t146 = cos(qJ(4));
t159 = t146 * t130 + t144 * t131;
t157 = qJD(1) * t145;
t156 = t147 * qJD(1);
t155 = qJD(1) * qJD(3);
t153 = -t144 * t130 + t146 * t131;
t152 = t147 * qJD(2) - t145 * t136;
t138 = -qJD(4) + t156;
t124 = t138 * qJ(5) - t159;
t151 = qJD(5) - t153;
t129 = -qJD(3) * pkin(3) - t152;
t135 = t144 * qJD(3) + t146 * t157;
t150 = -t135 * qJ(5) + t129;
t137 = t154 * qJD(1);
t134 = -t146 * qJD(3) + t144 * t157;
t125 = t134 * pkin(4) + t150;
t123 = t138 * pkin(4) + t151;
t122 = t160 * t134 + t150;
t121 = -t134 * pkin(5) + qJD(6) - t124;
t120 = t135 * pkin(5) + t160 * t138 + t151;
t1 = [t161, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t142 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t148, t145 ^ 2 * t161, t145 * t148 * t147, t145 * t155, t147 * t155, qJD(3) ^ 2 / 0.2e1, t152 * qJD(3) - t137 * t156, -t158 * qJD(3) + t137 * t157, t135 ^ 2 / 0.2e1, -t135 * t134, -t135 * t138, t134 * t138, t138 ^ 2 / 0.2e1, t129 * t134 - t153 * t138, t129 * t135 + t159 * t138, t123 * t135 + t124 * t134, -t123 * t138 - t125 * t134, t124 * t138 - t125 * t135, t125 ^ 2 / 0.2e1 + t124 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1, t120 * t135 - t121 * t134, -t121 * t138 - t122 * t135, t120 * t138 + t122 * t134, t122 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1 + t121 ^ 2 / 0.2e1;];
T_reg  = t1;
