% Calculate minimal parameter regressor of potential energy for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:56
% DurationCPUTime: 0.03s
% Computational Cost: add. (41->18), mult. (40->26), div. (0->0), fcn. (34->8), ass. (0->14)
t95 = g(3) * (qJ(2) + pkin(4));
t87 = qJ(1) + pkin(6);
t83 = sin(t87);
t85 = cos(t87);
t94 = g(1) * t85 + g(2) * t83;
t91 = sin(qJ(1));
t92 = cos(qJ(1));
t93 = -g(1) * t92 - g(2) * t91;
t89 = cos(pkin(7));
t88 = sin(pkin(7));
t86 = pkin(7) + qJ(4);
t84 = cos(t86);
t82 = sin(t86);
t1 = [0, t93, g(1) * t91 - g(2) * t92, t93 * pkin(1) - t95, -g(3) * t88 - t94 * t89, -g(3) * t89 + t94 * t88, -g(1) * t83 + g(2) * t85, -g(1) * (t92 * pkin(1) + t85 * pkin(2) + t83 * qJ(3)) - g(2) * (t91 * pkin(1) + t83 * pkin(2) - t85 * qJ(3)) - t95, 0, 0, 0, 0, 0, -g(3) * t82 - t94 * t84, -g(3) * t84 + t94 * t82;];
U_reg = t1;
