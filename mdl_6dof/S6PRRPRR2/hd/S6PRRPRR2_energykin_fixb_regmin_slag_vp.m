% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:46:36
% EndTime: 2021-01-16 03:46:36
% DurationCPUTime: 0.23s
% Computational Cost: add. (317->49), mult. (764->105), div. (0->0), fcn. (580->12), ass. (0->46)
t87 = qJD(2) ^ 2;
t100 = t87 / 0.2e1;
t99 = cos(qJ(5));
t83 = sin(qJ(2));
t96 = qJD(1) * sin(pkin(6));
t70 = qJD(2) * pkin(8) + t83 * t96;
t85 = cos(qJ(3));
t95 = qJD(1) * cos(pkin(6));
t74 = t85 * t95;
t82 = sin(qJ(3));
t92 = qJ(4) * qJD(2);
t59 = qJD(3) * pkin(3) + t74 + (-t70 - t92) * t82;
t97 = t85 * t70 + t82 * t95;
t60 = t85 * t92 + t97;
t76 = sin(pkin(12));
t78 = cos(pkin(12));
t50 = t76 * t59 + t78 * t60;
t48 = qJD(3) * pkin(9) + t50;
t86 = cos(qJ(2));
t90 = t86 * t96;
t65 = -t90 + qJD(4) + (-pkin(3) * t85 - pkin(2)) * qJD(2);
t93 = qJD(2) * t85;
t94 = qJD(2) * t82;
t67 = t76 * t94 - t78 * t93;
t68 = (t76 * t85 + t78 * t82) * qJD(2);
t55 = t67 * pkin(4) - t68 * pkin(9) + t65;
t81 = sin(qJ(5));
t98 = t99 * t48 + t81 * t55;
t91 = qJD(2) * qJD(3);
t89 = qJD(2) * t96;
t88 = -t81 * t48 + t99 * t55;
t49 = t78 * t59 - t76 * t60;
t66 = qJD(5) + t67;
t47 = -qJD(3) * pkin(4) - t49;
t84 = cos(qJ(6));
t80 = sin(qJ(6));
t71 = -qJD(2) * pkin(2) - t90;
t64 = qJD(6) + t66;
t63 = t81 * qJD(3) + t99 * t68;
t62 = -t99 * qJD(3) + t81 * t68;
t52 = -t80 * t62 + t84 * t63;
t51 = t84 * t62 + t80 * t63;
t45 = t62 * pkin(5) + t47;
t44 = -t62 * pkin(10) + t98;
t43 = t66 * pkin(5) - t63 * pkin(10) + t88;
t1 = [qJD(1) ^ 2 / 0.2e1, t100, t86 * t89, -t83 * t89, t82 ^ 2 * t100, t82 * t87 * t85, t82 * t91, t85 * t91, qJD(3) ^ 2 / 0.2e1, (-t82 * t70 + t74) * qJD(3) - t71 * t93, -t97 * qJD(3) + t71 * t94, t49 * qJD(3) + t65 * t67, -t50 * qJD(3) + t65 * t68, -t49 * t68 - t50 * t67, t50 ^ 2 / 0.2e1 + t49 ^ 2 / 0.2e1 + t65 ^ 2 / 0.2e1, t63 ^ 2 / 0.2e1, -t63 * t62, t63 * t66, -t62 * t66, t66 ^ 2 / 0.2e1, t47 * t62 + t88 * t66, t47 * t63 - t98 * t66, t52 ^ 2 / 0.2e1, -t52 * t51, t52 * t64, -t51 * t64, t64 ^ 2 / 0.2e1, (t84 * t43 - t80 * t44) * t64 + t45 * t51, -(t80 * t43 + t84 * t44) * t64 + t45 * t52;];
T_reg = t1;
