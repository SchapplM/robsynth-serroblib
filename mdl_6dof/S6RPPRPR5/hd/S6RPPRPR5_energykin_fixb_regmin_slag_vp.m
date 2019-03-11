% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:20
% EndTime: 2019-03-09 01:49:20
% DurationCPUTime: 0.07s
% Computational Cost: add. (176->39), mult. (324->84), div. (0->0), fcn. (163->6), ass. (0->32)
t136 = qJD(1) ^ 2;
t144 = t136 / 0.2e1;
t143 = -pkin(1) - qJ(3);
t142 = cos(pkin(9));
t133 = sin(qJ(4));
t135 = cos(qJ(4));
t116 = -qJD(2) + (pkin(4) * t133 - qJ(5) * t135 - t143) * qJD(1);
t127 = qJD(1) * qJ(2) + qJD(3);
t123 = -qJD(1) * pkin(7) + t127;
t120 = qJD(4) * qJ(5) + t133 * t123;
t131 = sin(pkin(9));
t110 = t131 * t116 + t142 * t120;
t141 = qJD(1) * t135;
t140 = qJD(4) * t123;
t124 = -t143 * qJD(1) - qJD(2);
t139 = t124 * qJD(1);
t138 = t133 * qJD(1);
t137 = qJD(1) * qJD(4);
t109 = t142 * t116 - t131 * t120;
t118 = -qJD(4) * pkin(4) - t135 * t123 + qJD(5);
t134 = cos(qJ(6));
t132 = sin(qJ(6));
t128 = -qJD(1) * pkin(1) + qJD(2);
t126 = qJD(6) + t138;
t122 = t131 * qJD(4) + t142 * t141;
t121 = -t142 * qJD(4) + t131 * t141;
t113 = t121 * pkin(5) + t118;
t112 = -t132 * t121 + t134 * t122;
t111 = t134 * t121 + t132 * t122;
t108 = -t121 * pkin(8) + t110;
t107 = pkin(5) * t138 - t122 * pkin(8) + t109;
t1 = [t144, 0, 0, t128 * qJD(1), t136 * qJ(2), qJ(2) ^ 2 * t144 + t128 ^ 2 / 0.2e1, t127 * qJD(1), t139, t124 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1, t135 ^ 2 * t144, -t135 * t136 * t133, t135 * t137, -t133 * t137, qJD(4) ^ 2 / 0.2e1, t124 * t138 + t135 * t140, -t133 * t140 + t135 * t139, t109 * t138 + t118 * t121, -t110 * t138 + t118 * t122, -t109 * t122 - t110 * t121, t110 ^ 2 / 0.2e1 + t109 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1, t112 ^ 2 / 0.2e1, -t112 * t111, t112 * t126, -t111 * t126, t126 ^ 2 / 0.2e1 (t134 * t107 - t132 * t108) * t126 + t113 * t111 -(t132 * t107 + t134 * t108) * t126 + t113 * t112;];
T_reg  = t1;
