% Calculate potential energy for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:01
% EndTime: 2019-12-31 17:49:02
% DurationCPUTime: 0.27s
% Computational Cost: add. (137->48), mult. (111->38), div. (0->0), fcn. (80->8), ass. (0->22)
t66 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t65 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t64 = -m(3) - m(4);
t63 = -m(5) - m(6);
t47 = pkin(8) + qJ(4);
t40 = sin(t47);
t42 = cos(t47);
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t62 = -m(4) * pkin(2) - t50 * mrSges(4,1) + t49 * mrSges(4,2) - t40 * t65 + t42 * t66 - mrSges(3,1);
t61 = m(4) * qJ(3) - mrSges(3,2) + mrSges(6,2) + mrSges(4,3) + mrSges(5,3);
t53 = sin(qJ(1));
t45 = t53 * pkin(1);
t54 = cos(qJ(1));
t46 = t54 * pkin(1);
t51 = qJ(2) + pkin(5);
t52 = -pkin(6) - qJ(3);
t48 = qJ(1) + pkin(7);
t43 = cos(t48);
t41 = sin(t48);
t39 = pkin(3) * t50 + pkin(2);
t1 = (-m(2) * pkin(5) - t49 * mrSges(4,1) - t50 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t63 * (t49 * pkin(3) + t51) + t64 * t51 + t65 * t42 + t66 * t40) * g(3) + (-t53 * mrSges(2,1) - t54 * mrSges(2,2) - mrSges(1,2) + t63 * (t41 * t39 + t43 * t52 + t45) + t64 * t45 + t61 * t43 + t62 * t41) * g(2) + (-t54 * mrSges(2,1) + t53 * mrSges(2,2) - mrSges(1,1) + t63 * (t43 * t39 - t41 * t52 + t46) + t64 * t46 + t62 * t43 - t61 * t41) * g(1);
U = t1;
