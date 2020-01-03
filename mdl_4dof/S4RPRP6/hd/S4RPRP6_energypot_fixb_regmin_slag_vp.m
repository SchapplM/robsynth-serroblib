% Calculate minimal parameter regressor of potential energy for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:13
% EndTime: 2019-12-31 16:46:13
% DurationCPUTime: 0.03s
% Computational Cost: add. (25->21), mult. (41->24), div. (0->0), fcn. (35->4), ass. (0->11)
t63 = sin(qJ(3));
t68 = pkin(3) * t63;
t64 = sin(qJ(1));
t66 = cos(qJ(1));
t67 = t66 * pkin(1) + t64 * qJ(2);
t57 = g(1) * t64 - g(2) * t66;
t65 = cos(qJ(3));
t62 = -qJ(4) - pkin(5);
t60 = t64 * pkin(1);
t58 = g(1) * t66 + g(2) * t64;
t1 = [0, -t58, t57, t58, -t57, -g(1) * t67 - g(2) * (-t66 * qJ(2) + t60) - g(3) * pkin(4), 0, 0, 0, 0, 0, -g(3) * t65 - t57 * t63, g(3) * t63 - t57 * t65, -t58, -g(1) * (-t66 * t62 + t64 * t68 + t67) - g(2) * (-t64 * t62 + t60 + (-qJ(2) - t68) * t66) - g(3) * (t65 * pkin(3) + pkin(2) + pkin(4));];
U_reg = t1;
