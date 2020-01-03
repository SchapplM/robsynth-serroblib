% Calculate potential energy for
% S5RPRRP5
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:38
% EndTime: 2019-12-31 18:40:38
% DurationCPUTime: 0.23s
% Computational Cost: add. (139->47), mult. (99->38), div. (0->0), fcn. (68->8), ass. (0->21)
t64 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t63 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t62 = -m(5) - m(6);
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t61 = -t63 * t49 + t64 * t51 - mrSges(4,1);
t60 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t50 = sin(qJ(1));
t45 = t50 * pkin(1);
t52 = cos(qJ(1));
t46 = t52 * pkin(1);
t59 = qJ(2) + pkin(5);
t48 = qJ(1) + pkin(8);
t42 = sin(t48);
t58 = pkin(2) * t42 + t45;
t43 = cos(t48);
t57 = pkin(2) * t43 + t46;
t44 = qJ(3) + t48;
t41 = cos(t44);
t40 = sin(t44);
t1 = (-m(2) * pkin(5) - m(3) * t59 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t63 * t51 + t64 * t49 + (-m(4) + t62) * (pkin(6) + t59)) * g(3) + (-m(3) * t45 - m(4) * t58 - t50 * mrSges(2,1) - t42 * mrSges(3,1) - t52 * mrSges(2,2) - t43 * mrSges(3,2) - mrSges(1,2) + t62 * (t40 * pkin(3) - t41 * pkin(7) + t58) - t60 * t41 + t61 * t40) * g(2) + (-m(3) * t46 - m(4) * t57 - t52 * mrSges(2,1) - t43 * mrSges(3,1) + t50 * mrSges(2,2) + t42 * mrSges(3,2) - mrSges(1,1) + t62 * (t41 * pkin(3) + t40 * pkin(7) + t57) + t61 * t41 + t60 * t40) * g(1);
U = t1;
