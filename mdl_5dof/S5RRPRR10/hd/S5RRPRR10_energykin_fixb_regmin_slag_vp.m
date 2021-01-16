% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:02:42
% EndTime: 2021-01-15 22:02:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (327->46), mult. (952->98), div. (0->0), fcn. (734->10), ass. (0->42)
t76 = sin(pkin(5));
t85 = qJD(1) ^ 2;
t95 = t76 ^ 2 * t85;
t84 = cos(qJ(2));
t91 = cos(pkin(5)) * qJD(1);
t89 = pkin(1) * t91;
t72 = t84 * t89;
t73 = qJD(2) + t91;
t81 = sin(qJ(2));
t92 = qJD(1) * t76;
t88 = t81 * t92;
t59 = t73 * pkin(2) + t72 + (-pkin(7) - qJ(3)) * t88;
t87 = t84 * t92;
t93 = pkin(7) * t87 + t81 * t89;
t62 = qJ(3) * t87 + t93;
t75 = sin(pkin(10));
t77 = cos(pkin(10));
t50 = t75 * t59 + t77 * t62;
t48 = t73 * pkin(8) + t50;
t65 = t75 * t88 - t77 * t87;
t66 = (t75 * t84 + t77 * t81) * t92;
t67 = qJD(3) + (-pkin(2) * t84 - pkin(1)) * t92;
t54 = t65 * pkin(3) - t66 * pkin(8) + t67;
t80 = sin(qJ(4));
t83 = cos(qJ(4));
t94 = t83 * t48 + t80 * t54;
t90 = t84 * t95;
t49 = t77 * t59 - t75 * t62;
t57 = t80 * t66 - t83 * t73;
t86 = -t80 * t48 + t83 * t54;
t47 = -t73 * pkin(3) - t49;
t82 = cos(qJ(5));
t79 = sin(qJ(5));
t64 = qJD(4) + t65;
t58 = t83 * t66 + t80 * t73;
t56 = qJD(5) + t57;
t52 = t82 * t58 + t79 * t64;
t51 = t79 * t58 - t82 * t64;
t45 = t57 * pkin(4) - t58 * pkin(9) + t47;
t44 = t64 * pkin(9) + t94;
t43 = -t64 * pkin(4) - t86;
t1 = [t85 / 0.2e1, 0, 0, t81 ^ 2 * t95 / 0.2e1, t81 * t90, t73 * t88, t73 * t87, t73 ^ 2 / 0.2e1, pkin(1) * t90 + (-pkin(7) * t88 + t72) * t73, -pkin(1) * t81 * t95 - t93 * t73, t49 * t73 + t67 * t65, -t50 * t73 + t67 * t66, -t49 * t66 - t50 * t65, t50 ^ 2 / 0.2e1 + t49 ^ 2 / 0.2e1 + t67 ^ 2 / 0.2e1, t58 ^ 2 / 0.2e1, -t58 * t57, t58 * t64, -t57 * t64, t64 ^ 2 / 0.2e1, t47 * t57 + t86 * t64, t47 * t58 - t94 * t64, t52 ^ 2 / 0.2e1, -t52 * t51, t52 * t56, -t51 * t56, t56 ^ 2 / 0.2e1, (-t79 * t44 + t82 * t45) * t56 + t43 * t51, -(t82 * t44 + t79 * t45) * t56 + t43 * t52;];
T_reg = t1;
