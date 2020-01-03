% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP1
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
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:11
% EndTime: 2020-01-03 11:59:11
% DurationCPUTime: 0.04s
% Computational Cost: add. (89->23), mult. (155->55), div. (0->0), fcn. (72->6), ass. (0->23)
t86 = qJD(1) + qJD(2);
t85 = t86 ^ 2;
t100 = t85 / 0.2e1;
t97 = pkin(1) * qJD(1);
t93 = cos(qJ(2)) * t97;
t80 = t86 * pkin(2) + t93;
t87 = sin(pkin(8));
t88 = cos(pkin(8));
t94 = sin(qJ(2)) * t97;
t77 = t88 * t80 - t87 * t94;
t99 = (-t86 * pkin(3) - t77) * t86;
t78 = t87 * t80 + t88 * t94;
t76 = t86 * pkin(7) + t78;
t89 = sin(qJ(4));
t91 = cos(qJ(4));
t98 = t89 * qJD(3) + t91 * t76;
t96 = qJ(5) * t86;
t95 = qJD(4) * t86;
t84 = t91 * qJD(3);
t73 = qJD(5) + (-pkin(4) * t91 - pkin(3)) * t86 - t77;
t72 = t91 * t96 + t98;
t71 = qJD(4) * pkin(4) + t84 + (-t76 - t96) * t89;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t100, t86 * t93, -t86 * t94, t78 ^ 2 / 0.2e1 + t77 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t89 ^ 2 * t100, t89 * t85 * t91, t89 * t95, t91 * t95, qJD(4) ^ 2 / 0.2e1, -t91 * t99 + (-t89 * t76 + t84) * qJD(4), -t98 * qJD(4) + t89 * t99, (-t71 * t89 + t72 * t91) * t86, t72 ^ 2 / 0.2e1 + t71 ^ 2 / 0.2e1 + t73 ^ 2 / 0.2e1;];
T_reg = t1;
