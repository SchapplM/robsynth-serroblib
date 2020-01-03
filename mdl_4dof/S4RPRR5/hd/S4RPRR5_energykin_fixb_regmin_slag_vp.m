% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:39
% EndTime: 2019-12-31 16:51:39
% DurationCPUTime: 0.03s
% Computational Cost: add. (42->13), mult. (72->36), div. (0->0), fcn. (21->4), ass. (0->19)
t61 = -qJD(1) + qJD(3);
t60 = t61 ^ 2;
t74 = t60 / 0.2e1;
t66 = qJD(1) ^ 2;
t73 = t66 / 0.2e1;
t57 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t63 = sin(qJ(3));
t65 = cos(qJ(3));
t68 = qJ(2) * qJD(1);
t67 = t65 * t57 - t63 * t68;
t72 = (-t61 * pkin(3) - t67) * t61;
t71 = t63 * t57 + t65 * t68;
t62 = sin(qJ(4));
t70 = qJD(4) * t62;
t64 = cos(qJ(4));
t69 = qJD(4) * t64;
t59 = -qJD(1) * pkin(1) + qJD(2);
t55 = t61 * pkin(6) + t71;
t1 = [t73, 0, 0, -t59 * qJD(1), t66 * qJ(2), qJ(2) ^ 2 * t73 + t59 ^ 2 / 0.2e1, t74, t67 * t61, -t71 * t61, t62 ^ 2 * t74, t62 * t60 * t64, t61 * t70, t61 * t69, qJD(4) ^ 2 / 0.2e1, -t55 * t70 - t64 * t72, -t55 * t69 + t62 * t72;];
T_reg = t1;
