% Calculate potential energy for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:19
% EndTime: 2019-12-05 18:01:19
% DurationCPUTime: 0.28s
% Computational Cost: add. (130->45), mult. (92->34), div. (0->0), fcn. (61->8), ass. (0->20)
t57 = m(5) + m(6);
t56 = -mrSges(5,2) - mrSges(6,2);
t55 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t54 = -m(4) - t57;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t53 = t57 * pkin(3) + t56 * t42 + t55 * t44 + mrSges(4,1);
t52 = -m(5) * pkin(7) + m(6) * (-qJ(5) - pkin(7)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t43 = sin(qJ(1));
t51 = t43 * pkin(1);
t45 = cos(qJ(1));
t38 = t45 * pkin(1);
t50 = qJ(2) + pkin(5);
t40 = qJ(1) + pkin(8);
t37 = qJ(3) + t40;
t36 = cos(t40);
t35 = sin(t40);
t33 = cos(t37);
t32 = sin(t37);
t1 = (-m(3) * t38 - t45 * mrSges(2,1) - t36 * mrSges(3,1) + t43 * mrSges(2,2) + t35 * mrSges(3,2) - mrSges(1,3) + t54 * (pkin(2) * t36 + t38) - t53 * t33 + t52 * t32) * g(3) + (m(3) * t51 + t43 * mrSges(2,1) + t35 * mrSges(3,1) + t45 * mrSges(2,2) + t36 * mrSges(3,2) - mrSges(1,2) + t54 * (-pkin(2) * t35 - t51) + t52 * t33 + t53 * t32) * g(2) + (-m(2) * pkin(5) - m(3) * t50 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t56 * t44 - t55 * t42 + t54 * (pkin(6) + t50)) * g(1);
U = t1;
