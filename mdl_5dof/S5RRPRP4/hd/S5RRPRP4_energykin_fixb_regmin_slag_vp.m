% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:55
% EndTime: 2019-12-31 19:52:55
% DurationCPUTime: 0.06s
% Computational Cost: add. (108->23), mult. (136->54), div. (0->0), fcn. (44->4), ass. (0->20)
t75 = qJD(1) + qJD(2);
t74 = t75 ^ 2;
t88 = t74 / 0.2e1;
t76 = sin(qJ(4));
t87 = t75 * t76;
t78 = cos(qJ(4));
t86 = t75 * t78;
t85 = pkin(1) * qJD(1);
t84 = qJD(4) * t76;
t83 = qJD(4) * t78;
t82 = sin(qJ(2)) * t85;
t81 = cos(qJ(2)) * t85;
t80 = qJD(3) - t81;
t72 = t75 * qJ(3) + t82;
t71 = -t75 * pkin(2) + t80;
t70 = (-pkin(2) - pkin(7)) * t75 + t80;
t69 = qJD(4) * qJ(5) + t76 * t70;
t68 = -qJD(4) * pkin(4) - t78 * t70 + qJD(5);
t67 = t82 + (pkin(4) * t76 - qJ(5) * t78 + qJ(3)) * t75;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t88, t75 * t81, -t75 * t82, t71 * t75, t72 * t75, t72 ^ 2 / 0.2e1 + t71 ^ 2 / 0.2e1, t78 ^ 2 * t88, -t78 * t74 * t76, t75 * t83, -t75 * t84, qJD(4) ^ 2 / 0.2e1, t70 * t83 + t72 * t87, -t70 * t84 + t72 * t86, -t68 * qJD(4) + t67 * t87, (t68 * t78 - t69 * t76) * t75, t69 * qJD(4) - t67 * t86, t69 ^ 2 / 0.2e1 + t67 ^ 2 / 0.2e1 + t68 ^ 2 / 0.2e1;];
T_reg = t1;
