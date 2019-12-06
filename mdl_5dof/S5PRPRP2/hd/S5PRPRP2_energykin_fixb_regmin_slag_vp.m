% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:03
% EndTime: 2019-12-05 15:31:03
% DurationCPUTime: 0.10s
% Computational Cost: add. (66->26), mult. (165->60), div. (0->0), fcn. (88->4), ass. (0->23)
t100 = cos(qJ(4));
t97 = sin(pkin(8));
t98 = cos(pkin(8));
t87 = qJD(3) + (-pkin(3) * t98 - pkin(6) * t97 - pkin(2)) * qJD(2);
t105 = qJ(3) * qJD(2);
t91 = t97 * qJD(1) + t98 * t105;
t99 = sin(qJ(4));
t109 = t100 * t91 + t99 * t87;
t101 = qJD(2) ^ 2;
t108 = t101 * t97 ^ 2;
t107 = qJD(2) * t97;
t106 = t98 * qJD(2);
t104 = t99 * t107;
t103 = t100 * t107;
t102 = t100 * t87 - t99 * t91;
t95 = t98 * qJD(1);
t94 = -qJD(2) * pkin(2) + qJD(3);
t92 = -qJD(4) + t106;
t89 = t97 * t105 - t95;
t84 = qJD(5) - t95 + (pkin(4) * t99 + qJ(3)) * t107;
t83 = -qJ(5) * t104 + t109;
t82 = -t92 * pkin(4) - qJ(5) * t103 + t102;
t1 = [qJD(1) ^ 2 / 0.2e1, t101 / 0.2e1, 0, 0, -t94 * t106, t94 * t107, (t89 * t97 + t91 * t98) * qJD(2), t91 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1 + t94 ^ 2 / 0.2e1, t100 ^ 2 * t108 / 0.2e1, -t100 * t99 * t108, -t92 * t103, t92 * t104, t92 ^ 2 / 0.2e1, -t102 * t92 + t89 * t104, t89 * t103 + t109 * t92, (-t100 * t82 - t83 * t99) * t107, t83 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1 + t84 ^ 2 / 0.2e1;];
T_reg = t1;
