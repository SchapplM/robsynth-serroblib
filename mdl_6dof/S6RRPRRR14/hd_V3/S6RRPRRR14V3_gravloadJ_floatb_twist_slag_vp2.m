% Calculate Gravitation load on the joints for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:27
% EndTime: 2019-04-12 15:03:29
% DurationCPUTime: 0.56s
% Computational Cost: add. (270->73), mult. (705->109), div. (0->0), fcn. (754->10), ass. (0->45)
t32 = cos(qJ(2));
t27 = sin(qJ(2));
t45 = m(4) + m(5) + m(6) + m(7);
t69 = t45 * qJ(3) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t36 = t69 * t27;
t50 = mrSges(3,1) + mrSges(4,1);
t76 = t50 * t32 + mrSges(2,1) + t36;
t48 = mrSges(6,2) - mrSges(7,3);
t24 = sin(qJ(6));
t29 = cos(qJ(6));
t40 = mrSges(7,1) * t29 - mrSges(7,2) * t24 + mrSges(6,1);
t28 = sin(qJ(1));
t33 = cos(qJ(1));
t72 = g(1) * t33 + g(2) * t28;
t49 = mrSges(5,2) - mrSges(6,3);
t67 = -mrSges(7,1) * t24 - mrSges(7,2) * t29 + t49;
t26 = sin(qJ(4));
t31 = cos(qJ(4));
t70 = mrSges(5,1) * t31 - t26 * t67 + t50;
t25 = sin(qJ(5));
t30 = cos(qJ(5));
t66 = -t48 * t25 + t40 * t30 + mrSges(5,1);
t60 = t25 * t27;
t59 = t26 * t27;
t58 = t27 * t30;
t57 = t27 * t31;
t56 = t27 * t33;
t55 = t28 * t26;
t53 = t31 * t32;
t52 = t31 * t33;
t51 = t33 * t26;
t47 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2);
t16 = t28 * t53 - t51;
t42 = -t16 * t25 + t28 * t58;
t4 = t16 * t30 + t28 * t60;
t14 = -t25 * t32 + t30 * t57;
t41 = t25 * t57 + t30 * t32;
t20 = t32 * t52 + t55;
t19 = -t28 * t31 + t32 * t51;
t15 = t32 * t55 + t52;
t8 = t20 * t30 + t25 * t56;
t7 = t20 * t25 - t30 * t56;
t2 = t19 * t24 + t29 * t8;
t1 = t19 * t29 - t24 * t8;
t3 = [(-t20 * mrSges(5,1) - t8 * mrSges(6,1) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t49 * t19 + t47 * t28 - t76 * t33 + t48 * t7) * g(2) + (t16 * mrSges(5,1) - t67 * t15 + t76 * t28 + t47 * t33 + t40 * t4 + t48 * t42) * g(1) (t48 * (t25 * t53 - t58) - t40 * (t30 * t53 + t60) - t70 * t32 - t36) * g(3) + (t14 * t40 + t70 * t27 - t32 * t69 - t48 * t41) * t72 (g(3) * t32 - t72 * t27) * t45 (t66 * t26 + t67 * t31) * g(3) * t27 + (t66 * t15 + t67 * t16) * g(2) + (t66 * t19 + t67 * t20) * g(1) (t48 * t14 + t40 * t41) * g(3) + (t48 * t4 - t40 * t42) * g(2) + (t40 * t7 + t48 * t8) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t15 * t29 - t24 * t4) * mrSges(7,1) + (-t15 * t24 - t29 * t4) * mrSges(7,2)) - g(3) * ((-t14 * t24 + t29 * t59) * mrSges(7,1) + (-t14 * t29 - t24 * t59) * mrSges(7,2))];
taug  = t3(:);
