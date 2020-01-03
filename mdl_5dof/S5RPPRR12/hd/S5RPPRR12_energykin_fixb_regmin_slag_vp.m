% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:21
% DurationCPUTime: 0.07s
% Computational Cost: add. (116->31), mult. (268->73), div. (0->0), fcn. (156->6), ass. (0->29)
t116 = qJD(1) ^ 2;
t121 = t116 / 0.2e1;
t113 = sin(qJ(4));
t115 = cos(qJ(4));
t110 = sin(pkin(8));
t102 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t118 = -pkin(6) * qJD(1) + t102;
t95 = t118 * t110;
t111 = cos(pkin(8));
t96 = t118 * t111;
t120 = t113 * t96 + t115 * t95;
t119 = qJD(1) * t110;
t104 = qJD(1) * qJ(2) + qJD(3);
t100 = pkin(3) * t119 + t104;
t117 = -t113 * t95 + t115 * t96;
t98 = (t110 * t115 + t111 * t113) * qJD(1);
t114 = cos(qJ(5));
t112 = sin(qJ(5));
t108 = t111 ^ 2;
t107 = t110 ^ 2;
t105 = -qJD(1) * pkin(1) + qJD(2);
t99 = (-t110 * t113 + t111 * t115) * qJD(1);
t97 = qJD(5) + t98;
t92 = t112 * qJD(4) + t114 * t99;
t91 = -t114 * qJD(4) + t112 * t99;
t90 = t98 * pkin(4) - t99 * pkin(7) + t100;
t89 = qJD(4) * pkin(7) + t120;
t88 = -qJD(4) * pkin(4) - t117;
t1 = [t121, 0, 0, t105 * qJD(1), t116 * qJ(2), qJ(2) ^ 2 * t121 + t105 ^ 2 / 0.2e1, t104 * t119, t104 * t111 * qJD(1), (-t107 - t108) * t102 * qJD(1), t104 ^ 2 / 0.2e1 + (t107 / 0.2e1 + t108 / 0.2e1) * t102 ^ 2, t99 ^ 2 / 0.2e1, -t99 * t98, t99 * qJD(4), -t98 * qJD(4), qJD(4) ^ 2 / 0.2e1, t117 * qJD(4) + t100 * t98, -t120 * qJD(4) + t100 * t99, t92 ^ 2 / 0.2e1, -t92 * t91, t92 * t97, -t91 * t97, t97 ^ 2 / 0.2e1, (-t112 * t89 + t114 * t90) * t97 + t88 * t91, -(t112 * t90 + t114 * t89) * t97 + t88 * t92;];
T_reg = t1;
