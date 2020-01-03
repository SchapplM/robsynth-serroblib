% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRP6
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
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:09
% EndTime: 2019-12-31 17:19:09
% DurationCPUTime: 0.05s
% Computational Cost: add. (63->23), mult. (171->56), div. (0->0), fcn. (91->4), ass. (0->24)
t92 = qJD(1) ^ 2;
t102 = t92 / 0.2e1;
t101 = cos(qJ(3));
t91 = cos(qJ(2));
t100 = t91 * t92;
t90 = sin(qJ(2));
t79 = (-pkin(2) * t91 - pkin(6) * t90 - pkin(1)) * qJD(1);
t97 = t91 * qJD(1);
t84 = pkin(5) * t97 + qJD(2) * pkin(6);
t89 = sin(qJ(3));
t99 = t101 * t84 + t89 * t79;
t98 = qJD(1) * t90;
t96 = qJD(1) * qJD(2);
t95 = t90 * t96;
t94 = t91 * t96;
t93 = t101 * t79 - t89 * t84;
t83 = -qJD(2) * pkin(2) + pkin(5) * t98;
t85 = -qJD(3) + t97;
t81 = t89 * qJD(2) + t101 * t98;
t80 = -t101 * qJD(2) + t89 * t98;
t76 = t80 * pkin(3) + qJD(4) + t83;
t75 = -t80 * qJ(4) + t99;
t74 = -t85 * pkin(3) - t81 * qJ(4) + t93;
t1 = [t102, 0, 0, t90 ^ 2 * t102, t90 * t100, t95, t94, qJD(2) ^ 2 / 0.2e1, pkin(1) * t100 - pkin(5) * t95, -t92 * pkin(1) * t90 - pkin(5) * t94, t81 ^ 2 / 0.2e1, -t81 * t80, -t81 * t85, t80 * t85, t85 ^ 2 / 0.2e1, t83 * t80 - t93 * t85, t83 * t81 + t99 * t85, -t74 * t81 - t75 * t80, t75 ^ 2 / 0.2e1 + t74 ^ 2 / 0.2e1 + t76 ^ 2 / 0.2e1;];
T_reg = t1;
