% Calculate Gravitation load on the joints for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:37
% EndTime: 2019-12-31 19:52:38
% DurationCPUTime: 0.33s
% Computational Cost: add. (214->56), mult. (207->60), div. (0->0), fcn. (157->6), ass. (0->29)
t64 = -mrSges(5,2) + mrSges(6,3);
t63 = -m(5) - m(6);
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t62 = -t23 * mrSges(6,1) + t64 * t25;
t61 = -mrSges(3,1) + mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t60 = mrSges(3,2) - mrSges(4,3) + t62;
t22 = qJ(1) + qJ(2);
t19 = sin(t22);
t20 = cos(t22);
t58 = -g(1) * t19 + g(2) * t20;
t47 = t20 * t23;
t56 = -mrSges(5,1) * t47 + t60 * t20 + (t63 * (-pkin(2) - pkin(7)) - t61) * t19;
t48 = t19 * t23;
t55 = -mrSges(5,1) * t48 + t60 * t19 + t61 * t20;
t24 = sin(qJ(1));
t53 = pkin(1) * t24;
t26 = cos(qJ(1));
t21 = t26 * pkin(1);
t50 = t20 * pkin(2) + t19 * qJ(3);
t43 = qJ(5) * t25;
t17 = t20 * pkin(7);
t40 = t17 + t50;
t39 = t21 + t50;
t10 = t20 * qJ(3);
t38 = -pkin(2) * t19 + t10;
t30 = pkin(4) * t47 - t20 * t43 + t10;
t29 = pkin(4) * t48 - t19 * t43 + t40;
t1 = [(-mrSges(2,1) * t26 + t24 * mrSges(2,2) - m(3) * t21 - m(4) * t39 - m(5) * (t17 + t39) - m(6) * (t21 + t29) + t55) * g(2) + (t24 * mrSges(2,1) + mrSges(2,2) * t26 + m(3) * t53 - m(4) * (t38 - t53) - m(5) * (t10 - t53) - m(6) * (t30 - t53) + t56) * g(1), (-m(4) * t50 - m(5) * t40 - m(6) * t29 + t55) * g(2) + (-m(4) * t38 - m(5) * t10 - m(6) * t30 + t56) * g(1), t58 * (m(4) - t63), (mrSges(5,1) * t23 - m(6) * (-pkin(4) * t23 + t43) - t62) * g(3) + ((m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1)) * t25 + (m(6) * qJ(5) + t64) * t23) * t58, (-g(3) * t23 - t58 * t25) * m(6)];
taug = t1(:);
