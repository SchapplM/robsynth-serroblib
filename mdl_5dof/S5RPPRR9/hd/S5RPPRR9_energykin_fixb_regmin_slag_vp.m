% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR9
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
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:47
% EndTime: 2019-12-31 18:02:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (93->29), mult. (185->68), div. (0->0), fcn. (81->6), ass. (0->27)
t105 = qJD(1) ^ 2;
t112 = t105 / 0.2e1;
t102 = sin(qJ(4));
t104 = cos(qJ(4));
t110 = qJ(2) * qJD(1);
t92 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t98 = sin(pkin(8));
t99 = cos(pkin(8));
t88 = t99 * t110 + t98 * t92;
t86 = -qJD(1) * pkin(6) + t88;
t111 = t102 * qJD(3) + t104 * t86;
t109 = qJD(1) * t102;
t108 = t104 * qJD(1);
t107 = qJD(1) * qJD(4);
t87 = -t98 * t110 + t99 * t92;
t85 = qJD(1) * pkin(3) - t87;
t106 = t104 * qJD(3) - t102 * t86;
t103 = cos(qJ(5));
t101 = sin(qJ(5));
t96 = -qJD(1) * pkin(1) + qJD(2);
t93 = qJD(5) + t108;
t90 = t101 * qJD(4) - t103 * t109;
t89 = t103 * qJD(4) + t101 * t109;
t83 = (pkin(4) * t104 + pkin(7) * t102) * qJD(1) + t85;
t82 = qJD(4) * pkin(7) + t111;
t81 = -qJD(4) * pkin(4) - t106;
t1 = [t112, 0, 0, -t96 * qJD(1), t105 * qJ(2), qJ(2) ^ 2 * t112 + t96 ^ 2 / 0.2e1, -t87 * qJD(1), t88 * qJD(1), t88 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t102 ^ 2 * t112, t102 * t105 * t104, -t102 * t107, -t104 * t107, qJD(4) ^ 2 / 0.2e1, t106 * qJD(4) + t85 * t108, -t111 * qJD(4) - t85 * t109, t90 ^ 2 / 0.2e1, t90 * t89, t90 * t93, t89 * t93, t93 ^ 2 / 0.2e1, (-t101 * t82 + t103 * t83) * t93 - t81 * t89, -(t101 * t83 + t103 * t82) * t93 + t81 * t90;];
T_reg = t1;
