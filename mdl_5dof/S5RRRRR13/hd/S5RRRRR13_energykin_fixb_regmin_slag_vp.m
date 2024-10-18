% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:31:09
% EndTime: 2024-09-27 17:31:09
% DurationCPUTime: 0.02s
% Computational Cost: add. (233->28), mult. (300->72), div. (0->0), fcn. (183->10), ass. (0->37)
t147 = pkin(1) * qJD(1);
t124 = qJD(1) + qJD(2);
t122 = qJD(3) + t124;
t139 = cos(qJ(2)) * t147;
t118 = t124 * pkin(2) + t139;
t129 = sin(qJ(3));
t133 = cos(qJ(3));
t140 = sin(qJ(2)) * t147;
t135 = t133 * t118 - t129 * t140;
t111 = t122 * pkin(3) + t135;
t126 = cos(pkin(5));
t146 = t111 * t126;
t121 = t122 ^ 2;
t125 = sin(pkin(5));
t123 = t125 ^ 2;
t145 = t121 * t123;
t144 = t122 * t125;
t132 = cos(qJ(4));
t143 = t122 * t132;
t141 = t129 * t118 + t133 * t140;
t110 = pkin(9) * t144 + t141;
t128 = sin(qJ(4));
t142 = t132 * t110 + t128 * t146;
t138 = t111 * t122 * t123;
t137 = t128 * t144;
t136 = t125 * t143;
t119 = t126 * t122 + qJD(4);
t131 = cos(qJ(5));
t127 = sin(qJ(5));
t117 = qJD(5) + t119;
t113 = (t127 * t132 + t128 * t131) * t144;
t112 = t127 * t137 - t131 * t136;
t108 = t132 * t146;
t106 = (-pkin(4) * t143 - t111) * t125;
t105 = pkin(10) * t136 + t142;
t104 = t119 * pkin(4) + t108 + (-pkin(10) * t144 - t110) * t128;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t124 ^ 2 / 0.2e1, t124 * t139, -t124 * t140, t121 / 0.2e1, t135 * t122, -t141 * t122, t128 ^ 2 * t145 / 0.2e1, t128 * t132 * t145, t119 * t137, t119 * t136, t119 ^ 2 / 0.2e1, (-t128 * t110 + t108) * t119 + t132 * t138, -t142 * t119 - t128 * t138, t113 ^ 2 / 0.2e1, -t113 * t112, t113 * t117, -t112 * t117, t117 ^ 2 / 0.2e1, (t131 * t104 - t127 * t105) * t117 + t106 * t112, -(t127 * t104 + t131 * t105) * t117 + t106 * t113;];
T_reg = t1;
