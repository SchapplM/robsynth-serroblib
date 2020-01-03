% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:25
% EndTime: 2019-12-31 21:49:25
% DurationCPUTime: 0.10s
% Computational Cost: add. (149->22), mult. (193->56), div. (0->0), fcn. (87->6), ass. (0->24)
t86 = qJD(1) + qJD(2);
t85 = qJD(3) + t86;
t84 = t85 ^ 2;
t102 = t84 / 0.2e1;
t87 = sin(qJ(4));
t101 = t85 * t87;
t90 = cos(qJ(4));
t100 = t85 * t90;
t98 = pkin(1) * qJD(1);
t94 = cos(qJ(2)) * t98;
t81 = t86 * pkin(2) + t94;
t88 = sin(qJ(3));
t91 = cos(qJ(3));
t95 = sin(qJ(2)) * t98;
t99 = t88 * t81 + t91 * t95;
t97 = qJD(4) * t87;
t96 = qJD(4) * t90;
t93 = t91 * t81 - t88 * t95;
t79 = t85 * pkin(8) + t99;
t78 = -t85 * pkin(3) - t93;
t77 = qJD(4) * qJ(5) + t90 * t79;
t76 = -qJD(4) * pkin(4) + t87 * t79 + qJD(5);
t75 = (-pkin(4) * t90 - qJ(5) * t87 - pkin(3)) * t85 - t93;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t86 ^ 2 / 0.2e1, t86 * t94, -t86 * t95, t102, t93 * t85, -t99 * t85, t87 ^ 2 * t102, t87 * t84 * t90, t85 * t97, t85 * t96, qJD(4) ^ 2 / 0.2e1, -t78 * t100 - t79 * t97, t78 * t101 - t79 * t96, -t76 * qJD(4) - t75 * t100, (t76 * t87 + t77 * t90) * t85, t77 * qJD(4) - t75 * t101, t77 ^ 2 / 0.2e1 + t75 ^ 2 / 0.2e1 + t76 ^ 2 / 0.2e1;];
T_reg = t1;
