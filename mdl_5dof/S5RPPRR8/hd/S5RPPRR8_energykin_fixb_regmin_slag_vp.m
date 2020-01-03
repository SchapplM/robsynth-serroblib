% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:13
% EndTime: 2019-12-31 18:01:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (89->22), mult. (150->53), div. (0->0), fcn. (55->6), ass. (0->24)
t78 = -qJD(1) + qJD(4);
t77 = t78 ^ 2;
t92 = t77 / 0.2e1;
t85 = qJD(1) ^ 2;
t91 = t85 / 0.2e1;
t75 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t80 = cos(pkin(8));
t74 = t80 * t75;
t79 = sin(pkin(8));
t71 = t74 + (-qJ(2) * t79 - pkin(3)) * qJD(1);
t87 = qJ(2) * qJD(1);
t73 = t75 * t79 + t80 * t87;
t82 = sin(qJ(4));
t84 = cos(qJ(4));
t86 = t71 * t84 - t73 * t82;
t90 = (-pkin(4) * t78 - t86) * t78;
t89 = t82 * t71 + t84 * t73;
t88 = qJD(5) * t78;
t83 = cos(qJ(5));
t81 = sin(qJ(5));
t76 = -pkin(1) * qJD(1) + qJD(2);
t72 = -t79 * t87 + t74;
t68 = pkin(7) * t78 + t89;
t1 = [t91, 0, 0, -t76 * qJD(1), t85 * qJ(2), qJ(2) ^ 2 * t91 + t76 ^ 2 / 0.2e1, -t72 * qJD(1), t73 * qJD(1), t73 ^ 2 / 0.2e1 + t72 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t92, t86 * t78, -t89 * t78, t81 ^ 2 * t92, t81 * t77 * t83, t81 * t88, t83 * t88, qJD(5) ^ 2 / 0.2e1, -t83 * t90 + (qJD(3) * t83 - t68 * t81) * qJD(5), t81 * t90 - (qJD(3) * t81 + t68 * t83) * qJD(5);];
T_reg = t1;
