% Calculate Gravitation load on the joints for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:13
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:13:12
% EndTime: 2018-11-23 15:13:13
% DurationCPUTime: 0.87s
% Computational Cost: add. (1073->108), mult. (1317->140), div. (0->0), fcn. (1280->14), ass. (0->67)
t113 = m(7) * pkin(5);
t109 = -m(5) - m(7);
t103 = m(6) - t109;
t108 = -mrSges(6,1) - mrSges(7,1);
t102 = t108 - t113;
t114 = -m(7) * (-qJ(6) - pkin(9)) + mrSges(7,3) + mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t112 = -mrSges(4,2) + mrSges(5,3);
t107 = -mrSges(6,2) - mrSges(7,2);
t110 = m(6) * pkin(9) + pkin(3) * t103 + t114;
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t49 = sin(qJ(5));
t89 = t49 * t50;
t99 = t112 * t50 + t89 * t113 + t114 * t53 + mrSges(3,1);
t47 = sin(pkin(10));
t51 = sin(qJ(2));
t80 = pkin(6) + qJ(2);
t71 = cos(t80) / 0.2e1;
t81 = pkin(6) - qJ(2);
t74 = cos(t81);
t58 = t74 / 0.2e1 + t71;
t83 = cos(pkin(10));
t26 = t47 * t51 - t58 * t83;
t85 = qJ(4) * t50;
t93 = t26 * t53;
t106 = -pkin(3) * t93 - t26 * t85;
t29 = t47 * t58 + t51 * t83;
t92 = t29 * t53;
t105 = -pkin(3) * t92 - t29 * t85;
t72 = sin(t80);
t69 = t72 / 0.2e1;
t73 = sin(t81);
t70 = t73 / 0.2e1;
t39 = t69 + t70;
t91 = t39 * t53;
t104 = pkin(3) * t91 + t39 * t85;
t52 = cos(qJ(5));
t98 = -t103 * qJ(4) + t102 * t49 + t107 * t52 - t112;
t97 = m(6) * (pkin(4) + pkin(8)) + m(7) * (pkin(5) * t52 + pkin(4)) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t54 = cos(qJ(2));
t90 = t47 * t54;
t88 = t50 * t52;
t84 = cos(pkin(6));
t82 = sin(pkin(6));
t23 = t26 * pkin(2);
t57 = t70 - t72 / 0.2e1;
t28 = -t57 * t83 + t90;
t79 = pkin(8) * t28 - t23;
t24 = t29 * pkin(2);
t75 = t83 * t54;
t31 = t47 * t57 + t75;
t78 = pkin(8) * t31 - t24;
t38 = t39 * pkin(2);
t40 = t71 - t74 / 0.2e1;
t77 = -pkin(8) * t40 + t38;
t76 = t47 * t82;
t63 = t83 * t82;
t61 = t69 - t73 / 0.2e1;
t33 = -t40 * t53 + t50 * t84;
t32 = -t40 * t50 - t53 * t84;
t30 = -t47 * t61 + t75;
t27 = t61 * t83 + t90;
t16 = t30 * t53 + t50 * t76;
t15 = t30 * t50 - t53 * t76;
t14 = t27 * t53 - t50 * t63;
t13 = t27 * t50 + t53 * t63;
t1 = [(-m(2) - m(3) - m(4) - t103) * g(3) (-m(4) * t77 - m(6) * (pkin(9) * t91 + t104 + t38) + t109 * (t77 + t104) + t108 * (t39 * t89 - t40 * t52) + t107 * (t39 * t88 + t40 * t49) + t97 * t40 - t99 * t39) * g(3) + (-m(4) * t79 - m(6) * (-pkin(9) * t93 + t106 - t23) + t109 * (t79 + t106) + t108 * (-t26 * t89 + t28 * t52) + t107 * (-t26 * t88 - t28 * t49) - t97 * t28 + t99 * t26) * g(2) + (-m(4) * t78 - m(6) * (-pkin(9) * t92 + t105 - t24) + t108 * (-t29 * t89 + t31 * t52) + t107 * (-t29 * t88 - t31 * t49) + t109 * (t78 + t105) - t97 * t31 + t99 * t29) * g(1) (t110 * t32 + t98 * t33) * g(3) + (t110 * t13 + t98 * t14) * g(2) + (t110 * t15 + t98 * t16) * g(1), t103 * (-g(1) * t15 - g(2) * t13 - g(3) * t32) (t107 * (-t32 * t49 + t39 * t52) + t102 * (t32 * t52 + t39 * t49)) * g(3) + (t107 * (-t13 * t49 - t26 * t52) + t102 * (t13 * t52 - t26 * t49)) * g(2) + (t107 * (-t15 * t49 - t29 * t52) + t102 * (t15 * t52 - t29 * t49)) * g(1) (-g(1) * t16 - g(2) * t14 - g(3) * t33) * m(7)];
taug  = t1(:);
