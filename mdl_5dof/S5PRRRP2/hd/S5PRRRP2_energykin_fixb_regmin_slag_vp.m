% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:56
% EndTime: 2019-12-05 16:41:56
% DurationCPUTime: 0.07s
% Computational Cost: add. (76->19), mult. (120->50), div. (0->0), fcn. (48->4), ass. (0->19)
t73 = qJD(2) + qJD(3);
t72 = t73 ^ 2;
t86 = t72 / 0.2e1;
t74 = sin(qJ(4));
t85 = t73 * t74;
t76 = cos(qJ(4));
t84 = t73 * t76;
t82 = pkin(2) * qJD(2);
t80 = sin(qJ(3)) * t82;
t69 = t73 * pkin(7) + t80;
t83 = t74 * qJD(1) + t76 * t69;
t81 = qJD(4) * t73;
t79 = cos(qJ(3)) * t82;
t78 = t76 * qJD(1) - t74 * t69;
t70 = -t73 * pkin(3) - t79;
t67 = qJD(4) * qJ(5) + t83;
t66 = -qJD(4) * pkin(4) + qJD(5) - t78;
t65 = -t79 + (-pkin(4) * t76 - qJ(5) * t74 - pkin(3)) * t73;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t86, t73 * t79, -t73 * t80, t74 ^ 2 * t86, t74 * t72 * t76, t74 * t81, t76 * t81, qJD(4) ^ 2 / 0.2e1, t78 * qJD(4) - t70 * t84, -t83 * qJD(4) + t70 * t85, -t66 * qJD(4) - t65 * t84, (t66 * t74 + t67 * t76) * t73, t67 * qJD(4) - t65 * t85, t67 ^ 2 / 0.2e1 + t65 ^ 2 / 0.2e1 + t66 ^ 2 / 0.2e1;];
T_reg = t1;
