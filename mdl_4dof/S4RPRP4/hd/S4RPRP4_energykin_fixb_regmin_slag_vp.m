% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:53
% EndTime: 2019-12-31 16:43:53
% DurationCPUTime: 0.04s
% Computational Cost: add. (45->18), mult. (124->51), div. (0->0), fcn. (48->4), ass. (0->17)
t68 = qJD(1) ^ 2;
t76 = t68 / 0.2e1;
t64 = sin(pkin(6));
t61 = (pkin(1) * t64 + pkin(5)) * qJD(1);
t66 = sin(qJ(3));
t67 = cos(qJ(3));
t75 = t66 * qJD(2) + t67 * t61;
t65 = cos(pkin(6));
t71 = -pkin(1) * t65 - pkin(2);
t59 = (-pkin(3) * t67 - qJ(4) * t66 + t71) * qJD(1);
t74 = qJD(1) * t59;
t73 = t71 * t68;
t72 = qJD(1) * qJD(3);
t70 = t67 * qJD(2) - t66 * t61;
t58 = qJD(3) * qJ(4) + t75;
t57 = -qJD(3) * pkin(3) + qJD(4) - t70;
t1 = [t76, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t64 ^ 2 / 0.2e1 + t65 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t68, t66 ^ 2 * t76, t66 * t68 * t67, t66 * t72, t67 * t72, qJD(3) ^ 2 / 0.2e1, t70 * qJD(3) - t67 * t73, -t75 * qJD(3) + t66 * t73, -t57 * qJD(3) - t67 * t74, (t57 * t66 + t58 * t67) * qJD(1), t58 * qJD(3) - t66 * t74, t58 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1 + t57 ^ 2 / 0.2e1;];
T_reg = t1;
