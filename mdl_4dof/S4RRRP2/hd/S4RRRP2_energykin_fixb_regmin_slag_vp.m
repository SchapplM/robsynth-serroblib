% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRP2
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
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:09
% EndTime: 2019-12-31 17:13:09
% DurationCPUTime: 0.03s
% Computational Cost: add. (48->14), mult. (84->41), div. (0->0), fcn. (31->4), ass. (0->17)
t60 = qJD(1) + qJD(2);
t59 = t60 ^ 2;
t72 = t59 / 0.2e1;
t70 = pkin(1) * qJD(1);
t66 = cos(qJ(2)) * t70;
t71 = (-t60 * pkin(2) - t66) * t60;
t61 = sin(qJ(3));
t69 = qJD(3) * t61;
t63 = cos(qJ(3));
t68 = qJD(3) * t63;
t67 = sin(qJ(2)) * t70;
t57 = t60 * pkin(6) + t67;
t65 = qJ(4) * t60 + t57;
t56 = -t66 + qJD(4) + (-pkin(3) * t63 - pkin(2)) * t60;
t55 = t65 * t63;
t54 = qJD(3) * pkin(3) - t65 * t61;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t72, t60 * t66, -t60 * t67, t61 ^ 2 * t72, t61 * t59 * t63, t60 * t69, t60 * t68, qJD(3) ^ 2 / 0.2e1, -t57 * t69 - t63 * t71, -t57 * t68 + t61 * t71, (-t54 * t61 + t55 * t63) * t60, t55 ^ 2 / 0.2e1 + t54 ^ 2 / 0.2e1 + t56 ^ 2 / 0.2e1;];
T_reg = t1;
