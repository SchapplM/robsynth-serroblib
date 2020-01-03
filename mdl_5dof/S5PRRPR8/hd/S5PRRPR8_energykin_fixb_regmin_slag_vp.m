% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:41
% EndTime: 2019-12-31 17:42:41
% DurationCPUTime: 0.04s
% Computational Cost: add. (71->17), mult. (133->49), div. (0->0), fcn. (81->8), ass. (0->23)
t83 = qJD(2) + qJD(3);
t82 = t83 ^ 2;
t97 = t82 / 0.2e1;
t91 = cos(qJ(2));
t81 = qJD(2) * pkin(2) + t91 * qJD(1);
t87 = sin(qJ(3));
t90 = cos(qJ(3));
t88 = sin(qJ(2));
t95 = qJD(1) * t88;
t92 = t90 * t81 - t87 * t95;
t77 = t83 * pkin(3) + t92;
t79 = t87 * t81 + t90 * t95;
t84 = sin(pkin(9));
t85 = cos(pkin(9));
t74 = t85 * t77 - t84 * t79;
t96 = (-t83 * pkin(4) - t74) * t83;
t75 = t84 * t77 + t85 * t79;
t94 = qJD(5) * t83;
t93 = qJD(1) * qJD(2);
t89 = cos(qJ(5));
t86 = sin(qJ(5));
t73 = t83 * pkin(7) + t75;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t91 * t93, -t88 * t93, t97, t92 * t83, -t79 * t83, t75 ^ 2 / 0.2e1 + t74 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t86 ^ 2 * t97, t86 * t82 * t89, t86 * t94, t89 * t94, qJD(5) ^ 2 / 0.2e1, (t89 * qJD(4) - t86 * t73) * qJD(5) - t89 * t96, -(t86 * qJD(4) + t89 * t73) * qJD(5) + t86 * t96;];
T_reg = t1;
