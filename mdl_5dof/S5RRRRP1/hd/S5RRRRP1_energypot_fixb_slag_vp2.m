% Calculate potential energy for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:41
% EndTime: 2019-12-05 18:44:41
% DurationCPUTime: 0.20s
% Computational Cost: add. (132->47), mult. (117->34), div. (0->0), fcn. (86->8), ass. (0->20)
t69 = -m(5) - m(6);
t68 = mrSges(5,2) + mrSges(6,2);
t67 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t49 = qJ(2) + qJ(3);
t45 = qJ(4) + t49;
t39 = sin(t45);
t40 = cos(t45);
t52 = cos(qJ(2));
t41 = t52 * pkin(2) + pkin(1);
t42 = sin(t49);
t43 = cos(t49);
t50 = sin(qJ(2));
t66 = -m(3) * pkin(1) - m(4) * t41 - t52 * mrSges(3,1) - t43 * mrSges(4,1) + t50 * mrSges(3,2) + t42 * mrSges(4,2) - mrSges(2,1) + t69 * (pkin(3) * t43 + t41) + t67 * t40 + t68 * t39;
t54 = -pkin(7) - pkin(6);
t48 = -pkin(8) + t54;
t65 = mrSges(2,2) + m(6) * (-qJ(5) + t48) - mrSges(6,3) + m(5) * t48 - mrSges(5,3) + m(4) * t54 - mrSges(4,3) - m(3) * pkin(6) - mrSges(3,3);
t64 = t50 * pkin(2) + pkin(5);
t53 = cos(qJ(1));
t51 = sin(qJ(1));
t1 = (-m(4) * t64 - t50 * mrSges(3,1) - t42 * mrSges(4,1) - t52 * mrSges(3,2) - t43 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + t69 * (pkin(3) * t42 + t64) - t68 * t40 + t67 * t39 + (-m(2) - m(3)) * pkin(5)) * g(3) + (t66 * t51 - t65 * t53 - mrSges(1,2)) * g(2) + (t65 * t51 + t66 * t53 - mrSges(1,1)) * g(1);
U = t1;
