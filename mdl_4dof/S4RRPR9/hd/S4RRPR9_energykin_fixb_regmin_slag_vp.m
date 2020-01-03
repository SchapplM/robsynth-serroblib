% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:03
% EndTime: 2019-12-31 17:10:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (105->29), mult. (272->68), div. (0->0), fcn. (163->6), ass. (0->28)
t113 = qJD(1) ^ 2;
t121 = t113 / 0.2e1;
t112 = cos(qJ(2));
t117 = t112 * qJD(1);
t103 = pkin(5) * t117 + qJD(2) * qJ(3);
t108 = sin(pkin(7));
t120 = cos(pkin(7));
t110 = sin(qJ(2));
t98 = (-pkin(2) * t112 - qJ(3) * t110 - pkin(1)) * qJD(1);
t93 = t120 * t103 + t108 * t98;
t119 = t112 * t113;
t118 = qJD(1) * t110;
t116 = qJD(1) * qJD(2);
t115 = t110 * t116;
t114 = t112 * t116;
t92 = -t108 * t103 + t120 * t98;
t102 = -qJD(2) * pkin(2) + pkin(5) * t118 + qJD(3);
t111 = cos(qJ(4));
t109 = sin(qJ(4));
t104 = -qJD(4) + t117;
t100 = t108 * qJD(2) + t120 * t118;
t99 = -t120 * qJD(2) + t108 * t118;
t94 = t99 * pkin(3) + t102;
t91 = t111 * t100 - t109 * t99;
t90 = t109 * t100 + t111 * t99;
t89 = -t99 * pkin(6) + t93;
t88 = -pkin(3) * t117 - t100 * pkin(6) + t92;
t1 = [t121, 0, 0, t110 ^ 2 * t121, t110 * t119, t115, t114, qJD(2) ^ 2 / 0.2e1, pkin(1) * t119 - pkin(5) * t115, -t113 * pkin(1) * t110 - pkin(5) * t114, t102 * t99 - t92 * t117, t102 * t100 + t93 * t117, -t92 * t100 - t93 * t99, t93 ^ 2 / 0.2e1 + t92 ^ 2 / 0.2e1 + t102 ^ 2 / 0.2e1, t91 ^ 2 / 0.2e1, -t91 * t90, -t91 * t104, t90 * t104, t104 ^ 2 / 0.2e1, -(-t109 * t89 + t111 * t88) * t104 + t94 * t90, (t109 * t88 + t111 * t89) * t104 + t94 * t91;];
T_reg = t1;
