% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:48
% EndTime: 2019-12-31 19:27:48
% DurationCPUTime: 0.07s
% Computational Cost: add. (105->22), mult. (135->50), div. (0->0), fcn. (49->6), ass. (0->20)
t76 = qJD(1) + qJD(2);
t75 = t76 ^ 2;
t90 = t75 / 0.2e1;
t88 = pkin(1) * qJD(1);
t85 = cos(qJ(2)) * t88;
t84 = qJD(3) - t85;
t71 = (-pkin(2) - pkin(3)) * t76 + t84;
t86 = sin(qJ(2)) * t88;
t74 = t76 * qJ(3) + t86;
t78 = sin(pkin(8));
t79 = cos(pkin(8));
t68 = t79 * t71 - t78 * t74;
t89 = (t76 * pkin(4) - t68) * t76;
t69 = t78 * t71 + t79 * t74;
t87 = qJD(5) * t76;
t82 = cos(qJ(5));
t80 = sin(qJ(5));
t73 = -t76 * pkin(2) + t84;
t67 = -t76 * pkin(7) + t69;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t90, t76 * t85, -t76 * t86, -t73 * t76, t74 * t76, t74 ^ 2 / 0.2e1 + t73 ^ 2 / 0.2e1, -t68 * t76, t69 * t76, t69 ^ 2 / 0.2e1 + t68 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t80 ^ 2 * t90, t80 * t75 * t82, -t80 * t87, -t82 * t87, qJD(5) ^ 2 / 0.2e1, t82 * t89 + (qJD(4) * t82 - t67 * t80) * qJD(5), -t80 * t89 - (qJD(4) * t80 + t67 * t82) * qJD(5);];
T_reg = t1;
