% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRR9
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
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:28
% EndTime: 2019-12-31 16:56:28
% DurationCPUTime: 0.04s
% Computational Cost: add. (43->21), mult. (107->52), div. (0->0), fcn. (47->4), ass. (0->19)
t73 = qJD(1) ^ 2;
t78 = t73 / 0.2e1;
t77 = t73 * qJ(2);
t72 = cos(qJ(3));
t76 = qJD(1) * t72;
t65 = qJD(2) + (-pkin(1) - pkin(5)) * qJD(1);
t75 = qJD(3) * t65;
t74 = qJD(1) * qJD(3);
t71 = cos(qJ(4));
t70 = sin(qJ(3));
t69 = sin(qJ(4));
t67 = -qJD(1) * pkin(1) + qJD(2);
t66 = t70 * qJD(1) + qJD(4);
t64 = t69 * qJD(3) + t71 * t76;
t63 = -t71 * qJD(3) + t69 * t76;
t62 = -qJD(3) * pkin(3) - t72 * t65;
t61 = qJD(3) * pkin(6) + t70 * t65;
t60 = (pkin(3) * t70 - pkin(6) * t72 + qJ(2)) * qJD(1);
t1 = [t78, 0, 0, t67 * qJD(1), t77, qJ(2) ^ 2 * t78 + t67 ^ 2 / 0.2e1, t72 ^ 2 * t78, -t72 * t73 * t70, t72 * t74, -t70 * t74, qJD(3) ^ 2 / 0.2e1, t70 * t77 + t72 * t75, -t70 * t75 + t72 * t77, t64 ^ 2 / 0.2e1, -t64 * t63, t64 * t66, -t63 * t66, t66 ^ 2 / 0.2e1, (t71 * t60 - t69 * t61) * t66 + t62 * t63, -(t69 * t60 + t71 * t61) * t66 + t62 * t64;];
T_reg = t1;
