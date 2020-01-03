% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:55
% EndTime: 2019-12-31 17:11:55
% DurationCPUTime: 0.05s
% Computational Cost: add. (56->27), mult. (153->62), div. (0->0), fcn. (69->4), ass. (0->26)
t91 = qJD(1) ^ 2;
t102 = t91 / 0.2e1;
t101 = -pkin(2) - pkin(6);
t90 = cos(qJ(2));
t100 = t90 * t91;
t99 = qJD(1) * t90;
t88 = sin(qJ(2));
t98 = t88 * qJD(1);
t97 = pkin(5) * t98 + qJD(3);
t96 = qJD(2) * qJ(3);
t95 = qJD(1) * qJD(2);
t94 = t88 * t95;
t93 = t90 * t95;
t92 = -qJ(3) * t88 - pkin(1);
t89 = cos(qJ(4));
t87 = sin(qJ(4));
t85 = qJD(4) + t98;
t84 = -pkin(5) * t99 - t96;
t83 = -qJD(2) * pkin(2) + t97;
t82 = t89 * qJD(2) - t87 * t99;
t81 = t87 * qJD(2) + t89 * t99;
t80 = (-pkin(2) * t90 + t92) * qJD(1);
t79 = t96 + (pkin(3) + pkin(5)) * t99;
t78 = pkin(3) * t98 + t101 * qJD(2) + t97;
t77 = (t101 * t90 + t92) * qJD(1);
t1 = [t102, 0, 0, t88 ^ 2 * t102, t88 * t100, t94, t93, qJD(2) ^ 2 / 0.2e1, pkin(1) * t100 - pkin(5) * t94, -t91 * pkin(1) * t88 - pkin(5) * t93, (t83 * t88 - t84 * t90) * qJD(1), t83 * qJD(2) + t80 * t99, -t84 * qJD(2) - t80 * t98, t80 ^ 2 / 0.2e1 + t84 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1, t82 ^ 2 / 0.2e1, -t82 * t81, t82 * t85, -t81 * t85, t85 ^ 2 / 0.2e1, (-t87 * t77 + t89 * t78) * t85 + t79 * t81, -(t89 * t77 + t87 * t78) * t85 + t79 * t82;];
T_reg = t1;
