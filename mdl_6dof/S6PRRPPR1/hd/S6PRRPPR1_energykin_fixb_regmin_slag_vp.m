% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:03
% EndTime: 2021-01-16 02:06:03
% DurationCPUTime: 0.12s
% Computational Cost: add. (381->50), mult. (903->106), div. (0->0), fcn. (672->12), ass. (0->45)
t85 = qJD(2) ^ 2;
t96 = t85 / 0.2e1;
t81 = sin(qJ(2));
t93 = qJD(1) * sin(pkin(6));
t68 = qJD(2) * pkin(8) + t81 * t93;
t83 = cos(qJ(3));
t92 = qJD(1) * cos(pkin(6));
t72 = t83 * t92;
t80 = sin(qJ(3));
t89 = qJ(4) * qJD(2);
t58 = qJD(3) * pkin(3) + t72 + (-t68 - t89) * t80;
t95 = t83 * t68 + t80 * t92;
t59 = t83 * t89 + t95;
t75 = sin(pkin(11));
t77 = cos(pkin(11));
t49 = t75 * t58 + t77 * t59;
t47 = qJD(3) * qJ(5) + t49;
t84 = cos(qJ(2));
t87 = t84 * t93;
t63 = -t87 + qJD(4) + (-pkin(3) * t83 - pkin(2)) * qJD(2);
t90 = qJD(2) * t83;
t91 = qJD(2) * t80;
t65 = t75 * t91 - t77 * t90;
t66 = (t75 * t83 + t77 * t80) * qJD(2);
t54 = t65 * pkin(4) - t66 * qJ(5) + t63;
t74 = sin(pkin(12));
t94 = cos(pkin(12));
t43 = t94 * t47 + t74 * t54;
t88 = qJD(2) * qJD(3);
t86 = qJD(2) * t93;
t42 = -t74 * t47 + t94 * t54;
t48 = t77 * t58 - t75 * t59;
t46 = -qJD(3) * pkin(4) + qJD(5) - t48;
t82 = cos(qJ(6));
t79 = sin(qJ(6));
t69 = -qJD(2) * pkin(2) - t87;
t64 = qJD(6) + t65;
t62 = t74 * qJD(3) + t94 * t66;
t61 = -t94 * qJD(3) + t74 * t66;
t51 = -t79 * t61 + t82 * t62;
t50 = t82 * t61 + t79 * t62;
t44 = t61 * pkin(5) + t46;
t41 = -t61 * pkin(9) + t43;
t40 = t65 * pkin(5) - t62 * pkin(9) + t42;
t1 = [qJD(1) ^ 2 / 0.2e1, t96, t84 * t86, -t81 * t86, t80 ^ 2 * t96, t80 * t85 * t83, t80 * t88, t83 * t88, qJD(3) ^ 2 / 0.2e1, (-t80 * t68 + t72) * qJD(3) - t69 * t90, -t95 * qJD(3) + t69 * t91, t48 * qJD(3) + t63 * t65, -t49 * qJD(3) + t63 * t66, -t48 * t66 - t49 * t65, t49 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1 + t63 ^ 2 / 0.2e1, t42 * t65 + t46 * t61, -t43 * t65 + t46 * t62, -t42 * t62 - t43 * t61, t43 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1, t51 ^ 2 / 0.2e1, -t51 * t50, t51 * t64, -t50 * t64, t64 ^ 2 / 0.2e1, (t82 * t40 - t79 * t41) * t64 + t44 * t50, -(t79 * t40 + t82 * t41) * t64 + t44 * t51;];
T_reg = t1;
