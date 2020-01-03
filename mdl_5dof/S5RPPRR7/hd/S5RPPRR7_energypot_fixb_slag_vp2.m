% Calculate potential energy for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:28
% EndTime: 2019-12-31 17:59:28
% DurationCPUTime: 0.26s
% Computational Cost: add. (123->48), mult. (114->38), div. (0->0), fcn. (87->8), ass. (0->20)
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t66 = -m(6) * pkin(4) - t50 * mrSges(6,1) + t47 * mrSges(6,2) - mrSges(5,1);
t65 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t64 = -m(5) - m(6);
t62 = -t47 * mrSges(6,1) - t50 * mrSges(6,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t61 = t66 * t48 - t65 * t51 + mrSges(3,2) - mrSges(4,3);
t49 = sin(qJ(1));
t43 = t49 * pkin(1);
t52 = cos(qJ(1));
t44 = t52 * pkin(1);
t46 = qJ(2) + pkin(5);
t45 = qJ(1) + pkin(8);
t41 = sin(t45);
t60 = t41 * pkin(2) + t43;
t42 = cos(t45);
t57 = t42 * pkin(2) + t41 * qJ(3) + t44;
t1 = (-m(2) * pkin(5) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t64 * (pkin(3) + t46) + t66 * t51 + t65 * t48 + (-m(3) - m(4)) * t46) * g(3) + (-m(3) * t43 - m(4) * t60 - t49 * mrSges(2,1) - t52 * mrSges(2,2) - mrSges(1,2) + t64 * (t41 * pkin(6) + t60) + ((m(4) - t64) * qJ(3) - t61) * t42 + t62 * t41) * g(2) + (-m(3) * t44 - m(4) * t57 - t52 * mrSges(2,1) + t49 * mrSges(2,2) - mrSges(1,1) + t64 * (t42 * pkin(6) + t57) + t62 * t42 + t61 * t41) * g(1);
U = t1;
