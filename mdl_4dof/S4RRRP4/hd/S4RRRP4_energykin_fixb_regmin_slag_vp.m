% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRP4
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:41
% EndTime: 2019-12-31 17:15:41
% DurationCPUTime: 0.07s
% Computational Cost: add. (63->22), mult. (184->55), div. (0->0), fcn. (104->4), ass. (0->25)
t82 = qJD(1) ^ 2;
t93 = t82 / 0.2e1;
t92 = -pkin(6) - pkin(5);
t91 = cos(qJ(3));
t81 = cos(qJ(2));
t90 = t81 * t82;
t80 = sin(qJ(2));
t88 = qJD(1) * t80;
t74 = qJD(2) * pkin(2) + t92 * t88;
t87 = qJD(1) * t81;
t75 = t92 * t87;
t79 = sin(qJ(3));
t89 = t79 * t74 - t91 * t75;
t86 = qJD(1) * qJD(2);
t85 = t80 * t86;
t84 = t81 * t86;
t83 = t91 * t74 + t75 * t79;
t76 = (-pkin(2) * t81 - pkin(1)) * qJD(1);
t78 = qJD(2) + qJD(3);
t72 = (t79 * t81 + t91 * t80) * qJD(1);
t71 = t79 * t88 - t91 * t87;
t68 = pkin(3) * t71 + qJD(4) + t76;
t67 = -qJ(4) * t71 + t89;
t66 = pkin(3) * t78 - qJ(4) * t72 + t83;
t1 = [t93, 0, 0, t80 ^ 2 * t93, t80 * t90, t85, t84, qJD(2) ^ 2 / 0.2e1, pkin(1) * t90 - pkin(5) * t85, -pkin(1) * t80 * t82 - pkin(5) * t84, t72 ^ 2 / 0.2e1, -t72 * t71, t72 * t78, -t71 * t78, t78 ^ 2 / 0.2e1, t76 * t71 + t83 * t78, t76 * t72 - t89 * t78, -t66 * t72 - t67 * t71, t67 ^ 2 / 0.2e1 + t66 ^ 2 / 0.2e1 + t68 ^ 2 / 0.2e1;];
T_reg = t1;
