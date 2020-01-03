% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:41
% EndTime: 2019-12-31 16:52:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (86->29), mult. (259->65), div. (0->0), fcn. (172->6), ass. (0->27)
t105 = cos(qJ(3));
t104 = pkin(5) + qJ(2);
t93 = sin(pkin(7));
t102 = qJD(1) * t93;
t85 = t104 * t102;
t94 = cos(pkin(7));
t101 = qJD(1) * t94;
t86 = t104 * t101;
t96 = sin(qJ(3));
t103 = t105 * t86 - t96 * t85;
t100 = -t105 * t85 - t96 * t86;
t87 = qJD(2) + (-pkin(2) * t94 - pkin(1)) * qJD(1);
t98 = qJD(1) ^ 2;
t97 = cos(qJ(4));
t95 = sin(qJ(4));
t92 = qJD(3) + qJD(4);
t91 = t94 ^ 2;
t90 = t93 ^ 2;
t89 = -qJD(1) * pkin(1) + qJD(2);
t84 = (t105 * t93 + t94 * t96) * qJD(1);
t83 = -t105 * t101 + t96 * t102;
t78 = t83 * pkin(3) + t87;
t77 = -t95 * t83 + t97 * t84;
t76 = t97 * t83 + t95 * t84;
t75 = -t83 * pkin(6) + t103;
t74 = qJD(3) * pkin(3) - t84 * pkin(6) + t100;
t1 = [t98 / 0.2e1, 0, 0, -t89 * t101, t89 * t102, (t90 + t91) * t98 * qJ(2), t89 ^ 2 / 0.2e1 + (t91 / 0.2e1 + t90 / 0.2e1) * qJ(2) ^ 2 * t98, t84 ^ 2 / 0.2e1, -t84 * t83, t84 * qJD(3), -t83 * qJD(3), qJD(3) ^ 2 / 0.2e1, t100 * qJD(3) + t87 * t83, -t103 * qJD(3) + t87 * t84, t77 ^ 2 / 0.2e1, -t77 * t76, t77 * t92, -t76 * t92, t92 ^ 2 / 0.2e1, t78 * t76 + (t97 * t74 - t95 * t75) * t92, t78 * t77 - (t95 * t74 + t97 * t75) * t92;];
T_reg = t1;
