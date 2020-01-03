% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:18
% EndTime: 2019-12-31 17:23:18
% DurationCPUTime: 0.04s
% Computational Cost: add. (76->19), mult. (122->53), div. (0->0), fcn. (63->6), ass. (0->24)
t73 = qJD(1) + qJD(2);
t71 = t73 ^ 2;
t88 = t71 / 0.2e1;
t75 = sin(qJ(3));
t87 = t73 * t75;
t78 = cos(qJ(3));
t86 = t73 * t78;
t85 = pkin(1) * qJD(1);
t84 = qJD(3) * t75;
t83 = qJD(3) * t78;
t82 = sin(qJ(2)) * t85;
t81 = cos(qJ(2)) * t85;
t68 = t73 * pkin(6) + t82;
t80 = pkin(7) * t73 + t68;
t77 = cos(qJ(4));
t74 = sin(qJ(4));
t72 = qJD(3) + qJD(4);
t69 = -t73 * pkin(2) - t81;
t67 = -t81 + (-pkin(3) * t78 - pkin(2)) * t73;
t66 = (t74 * t78 + t75 * t77) * t73;
t65 = t74 * t87 - t77 * t86;
t64 = t80 * t78;
t63 = qJD(3) * pkin(3) - t80 * t75;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t88, t73 * t81, -t73 * t82, t75 ^ 2 * t88, t75 * t71 * t78, t73 * t84, t73 * t83, qJD(3) ^ 2 / 0.2e1, -t68 * t84 - t69 * t86, -t68 * t83 + t69 * t87, t66 ^ 2 / 0.2e1, -t66 * t65, t66 * t72, -t65 * t72, t72 ^ 2 / 0.2e1, t67 * t65 + (t77 * t63 - t74 * t64) * t72, t67 * t66 - (t74 * t63 + t77 * t64) * t72;];
T_reg = t1;
