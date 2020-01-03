% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:13
% EndTime: 2019-12-31 19:51:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (199->30), mult. (291->66), div. (0->0), fcn. (157->6), ass. (0->27)
t97 = qJD(1) + qJD(2);
t98 = sin(pkin(8));
t112 = t97 * t98;
t99 = cos(pkin(8));
t111 = t97 * t99;
t100 = sin(qJ(4));
t102 = cos(qJ(4));
t109 = pkin(1) * qJD(1);
t108 = sin(qJ(2)) * t109;
t93 = t97 * qJ(3) + t108;
t106 = pkin(7) * t97 + t93;
t86 = t106 * t98;
t87 = t106 * t99;
t110 = -t100 * t86 + t102 * t87;
t107 = cos(qJ(2)) * t109;
t105 = qJD(3) - t107;
t104 = -t100 * t87 - t102 * t86;
t88 = (-pkin(3) * t99 - pkin(2)) * t97 + t105;
t96 = t99 ^ 2;
t95 = t98 ^ 2;
t91 = -t97 * pkin(2) + t105;
t90 = (t100 * t99 + t102 * t98) * t97;
t89 = t100 * t112 - t102 * t111;
t83 = qJD(4) * qJ(5) + t110;
t82 = -qJD(4) * pkin(4) + qJD(5) - t104;
t81 = t89 * pkin(4) - t90 * qJ(5) + t88;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t97 ^ 2 / 0.2e1, t97 * t107, -t97 * t108, -t91 * t111, t91 * t112, (t95 + t96) * t97 * t93, t91 ^ 2 / 0.2e1 + (t96 / 0.2e1 + t95 / 0.2e1) * t93 ^ 2, t90 ^ 2 / 0.2e1, -t90 * t89, t90 * qJD(4), -t89 * qJD(4), qJD(4) ^ 2 / 0.2e1, t104 * qJD(4) + t88 * t89, -t110 * qJD(4) + t88 * t90, -t82 * qJD(4) + t81 * t89, t82 * t90 - t83 * t89, t83 * qJD(4) - t81 * t90, t83 ^ 2 / 0.2e1 + t81 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1;];
T_reg = t1;
