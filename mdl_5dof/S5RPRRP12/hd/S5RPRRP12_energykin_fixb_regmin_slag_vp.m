% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:19
% EndTime: 2019-12-31 18:57:19
% DurationCPUTime: 0.05s
% Computational Cost: add. (90->29), mult. (193->64), div. (0->0), fcn. (91->4), ass. (0->24)
t99 = qJD(1) ^ 2;
t107 = t99 / 0.2e1;
t106 = cos(qJ(4));
t97 = sin(qJ(3));
t98 = cos(qJ(3));
t87 = (pkin(3) * t97 - pkin(7) * t98 + qJ(2)) * qJD(1);
t92 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t88 = qJD(3) * pkin(7) + t97 * t92;
t96 = sin(qJ(4));
t105 = t106 * t88 + t96 * t87;
t104 = t99 * qJ(2);
t103 = qJD(1) * t98;
t102 = qJD(3) * t92;
t101 = qJD(1) * qJD(3);
t100 = t106 * t87 - t96 * t88;
t89 = -qJD(3) * pkin(3) - t98 * t92;
t94 = -qJD(1) * pkin(1) + qJD(2);
t93 = t97 * qJD(1) + qJD(4);
t91 = t96 * qJD(3) + t106 * t103;
t90 = -t106 * qJD(3) + t96 * t103;
t83 = t90 * pkin(4) + qJD(5) + t89;
t82 = -t90 * qJ(5) + t105;
t81 = t93 * pkin(4) - t91 * qJ(5) + t100;
t1 = [t107, 0, 0, t94 * qJD(1), t104, qJ(2) ^ 2 * t107 + t94 ^ 2 / 0.2e1, t98 ^ 2 * t107, -t98 * t99 * t97, t98 * t101, -t97 * t101, qJD(3) ^ 2 / 0.2e1, t98 * t102 + t97 * t104, -t97 * t102 + t98 * t104, t91 ^ 2 / 0.2e1, -t91 * t90, t91 * t93, -t90 * t93, t93 ^ 2 / 0.2e1, t100 * t93 + t89 * t90, -t105 * t93 + t89 * t91, -t81 * t91 - t82 * t90, t82 ^ 2 / 0.2e1 + t81 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1;];
T_reg = t1;
