% Calculate potential energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:34
% EndTime: 2019-12-05 17:41:34
% DurationCPUTime: 0.25s
% Computational Cost: add. (135->51), mult. (104->38), div. (0->0), fcn. (73->10), ass. (0->22)
t63 = -m(3) - m(4);
t64 = pkin(1) * (m(5) + m(6) - t63) + mrSges(2,1);
t47 = pkin(9) + qJ(4);
t42 = qJ(5) + t47;
t35 = sin(t42);
t36 = cos(t42);
t50 = cos(pkin(9));
t37 = t50 * pkin(3) + pkin(2);
t38 = sin(t47);
t40 = cos(t47);
t49 = sin(pkin(9));
t61 = mrSges(3,1) + m(4) * pkin(2) + t50 * mrSges(4,1) - t49 * mrSges(4,2) + m(5) * t37 + t40 * mrSges(5,1) - t38 * mrSges(5,2) + m(6) * (pkin(4) * t40 + t37) + t36 * mrSges(6,1) - t35 * mrSges(6,2);
t52 = -pkin(6) - qJ(3);
t60 = -m(4) * qJ(3) + m(5) * t52 + m(6) * (-pkin(7) + t52) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t51 = qJ(2) + pkin(5);
t58 = t49 * pkin(3) + t51;
t54 = cos(qJ(1));
t53 = sin(qJ(1));
t48 = qJ(1) + pkin(8);
t41 = cos(t48);
t39 = sin(t48);
t1 = (t53 * mrSges(2,2) + t60 * t39 - t61 * t41 - t64 * t54 - mrSges(1,3)) * g(3) + (t54 * mrSges(2,2) + t61 * t39 + t60 * t41 + t64 * t53 - mrSges(1,2)) * g(2) + (-mrSges(1,1) - m(2) * pkin(5) - mrSges(2,3) - mrSges(3,3) - t49 * mrSges(4,1) - t50 * mrSges(4,2) - m(5) * t58 - t38 * mrSges(5,1) - t40 * mrSges(5,2) - m(6) * (pkin(4) * t38 + t58) - t35 * mrSges(6,1) - t36 * mrSges(6,2) + t63 * t51) * g(1);
U = t1;
