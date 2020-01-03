% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:03
% EndTime: 2019-12-31 17:21:03
% DurationCPUTime: 0.05s
% Computational Cost: add. (91->25), mult. (221->60), div. (0->0), fcn. (118->4), ass. (0->24)
t100 = qJD(1) ^ 2;
t109 = t100 / 0.2e1;
t97 = sin(qJ(2));
t99 = cos(qJ(2));
t86 = (-pkin(2) * t99 - pkin(6) * t97 - pkin(1)) * qJD(1);
t105 = t99 * qJD(1);
t91 = pkin(5) * t105 + qJD(2) * pkin(6);
t96 = sin(qJ(3));
t98 = cos(qJ(3));
t108 = t96 * t86 + t98 * t91;
t107 = t100 * t99;
t106 = qJD(1) * t97;
t104 = qJD(1) * qJD(2);
t103 = t97 * t104;
t102 = t99 * t104;
t90 = -qJD(2) * pkin(2) + pkin(5) * t106;
t101 = t98 * t86 - t96 * t91;
t92 = -qJD(3) + t105;
t88 = t96 * qJD(2) + t98 * t106;
t87 = -t98 * qJD(2) + t96 * t106;
t84 = t87 * pkin(3) - t88 * qJ(4) + t90;
t83 = -t92 * qJ(4) + t108;
t82 = t92 * pkin(3) + qJD(4) - t101;
t1 = [t109, 0, 0, t97 ^ 2 * t109, t97 * t107, t103, t102, qJD(2) ^ 2 / 0.2e1, pkin(1) * t107 - pkin(5) * t103, -t100 * pkin(1) * t97 - pkin(5) * t102, t88 ^ 2 / 0.2e1, -t88 * t87, -t88 * t92, t87 * t92, t92 ^ 2 / 0.2e1, -t101 * t92 + t90 * t87, t108 * t92 + t90 * t88, t82 * t92 + t84 * t87, t82 * t88 - t83 * t87, -t83 * t92 - t84 * t88, t83 ^ 2 / 0.2e1 + t84 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1;];
T_reg = t1;
