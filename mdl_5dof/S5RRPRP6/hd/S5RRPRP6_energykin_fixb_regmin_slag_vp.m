% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:27
% EndTime: 2019-12-31 19:58:27
% DurationCPUTime: 0.07s
% Computational Cost: add. (167->34), mult. (426->75), div. (0->0), fcn. (274->6), ass. (0->34)
t140 = qJD(1) ^ 2;
t151 = t140 / 0.2e1;
t150 = cos(qJ(4));
t149 = pkin(6) + qJ(3);
t139 = cos(qJ(2));
t148 = t139 * t140;
t135 = sin(pkin(8));
t136 = cos(pkin(8));
t145 = qJD(1) * t139;
t138 = sin(qJ(2));
t146 = qJD(1) * t138;
t126 = -t135 * t146 + t136 * t145;
t127 = (t135 * t139 + t136 * t138) * qJD(1);
t132 = qJD(3) + (-pkin(2) * t139 - pkin(1)) * qJD(1);
t116 = -t126 * pkin(3) - t127 * pkin(7) + t132;
t130 = qJD(2) * pkin(2) - t149 * t146;
t131 = t149 * t145;
t121 = t135 * t130 + t136 * t131;
t119 = qJD(2) * pkin(7) + t121;
t137 = sin(qJ(4));
t147 = t137 * t116 + t150 * t119;
t144 = qJD(1) * qJD(2);
t143 = t138 * t144;
t142 = t139 * t144;
t141 = t150 * t116 - t137 * t119;
t120 = t136 * t130 - t135 * t131;
t118 = -qJD(2) * pkin(3) - t120;
t125 = qJD(4) - t126;
t123 = t137 * qJD(2) + t150 * t127;
t122 = -t150 * qJD(2) + t137 * t127;
t113 = t122 * pkin(4) + qJD(5) + t118;
t112 = -t122 * qJ(5) + t147;
t111 = t125 * pkin(4) - t123 * qJ(5) + t141;
t1 = [t151, 0, 0, t138 ^ 2 * t151, t138 * t148, t143, t142, qJD(2) ^ 2 / 0.2e1, pkin(1) * t148 - pkin(6) * t143, -t140 * pkin(1) * t138 - pkin(6) * t142, -t120 * t127 + t121 * t126, t121 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1, t123 ^ 2 / 0.2e1, -t123 * t122, t123 * t125, -t122 * t125, t125 ^ 2 / 0.2e1, t118 * t122 + t141 * t125, t118 * t123 - t147 * t125, -t111 * t123 - t112 * t122, t112 ^ 2 / 0.2e1 + t111 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1;];
T_reg = t1;
