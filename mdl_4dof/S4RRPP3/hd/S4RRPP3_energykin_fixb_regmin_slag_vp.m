% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:50
% EndTime: 2019-12-31 16:57:50
% DurationCPUTime: 0.05s
% Computational Cost: add. (85->23), mult. (229->56), div. (0->0), fcn. (121->4), ass. (0->24)
t91 = qJD(1) ^ 2;
t99 = t91 / 0.2e1;
t90 = cos(qJ(2));
t98 = t90 * t91;
t97 = pkin(5) + qJ(3);
t89 = sin(qJ(2));
t96 = qJD(1) * t89;
t83 = qJD(2) * pkin(2) - t97 * t96;
t95 = qJD(1) * t90;
t84 = t97 * t95;
t87 = sin(pkin(6));
t88 = cos(pkin(6));
t78 = t87 * t83 + t88 * t84;
t94 = qJD(1) * qJD(2);
t93 = t89 * t94;
t92 = t90 * t94;
t77 = t88 * t83 - t87 * t84;
t85 = qJD(3) + (-pkin(2) * t90 - pkin(1)) * qJD(1);
t81 = (t87 * t90 + t88 * t89) * qJD(1);
t80 = t87 * t96 - t88 * t95;
t76 = qJD(2) * qJ(4) + t78;
t75 = -qJD(2) * pkin(3) + qJD(4) - t77;
t74 = t80 * pkin(3) - t81 * qJ(4) + t85;
t1 = [t99, 0, 0, t89 ^ 2 * t99, t89 * t98, t93, t92, qJD(2) ^ 2 / 0.2e1, pkin(1) * t98 - pkin(5) * t93, -t91 * pkin(1) * t89 - pkin(5) * t92, -t77 * t81 - t78 * t80, t78 ^ 2 / 0.2e1 + t77 ^ 2 / 0.2e1 + t85 ^ 2 / 0.2e1, -t75 * qJD(2) + t74 * t80, t75 * t81 - t76 * t80, t76 * qJD(2) - t74 * t81, t76 ^ 2 / 0.2e1 + t74 ^ 2 / 0.2e1 + t75 ^ 2 / 0.2e1;];
T_reg = t1;
