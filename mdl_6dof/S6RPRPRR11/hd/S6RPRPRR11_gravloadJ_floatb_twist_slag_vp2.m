% Calculate Gravitation load on the joints for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:11:23
% EndTime: 2019-03-09 04:11:26
% DurationCPUTime: 1.20s
% Computational Cost: add. (1015->129), mult. (2428->186), div. (0->0), fcn. (3041->16), ass. (0->64)
t105 = cos(pkin(7));
t64 = cos(qJ(1));
t102 = sin(pkin(12));
t112 = sin(qJ(1));
t104 = cos(pkin(12));
t106 = cos(pkin(6));
t89 = t106 * t104;
t79 = t112 * t102 - t64 * t89;
t103 = sin(pkin(7));
t58 = sin(pkin(6));
t97 = t58 * t103;
t128 = t79 * t105 + t64 * t97;
t113 = cos(qJ(3));
t87 = t106 * t102;
t41 = t104 * t112 + t64 * t87;
t62 = sin(qJ(3));
t22 = -t113 * t41 + t128 * t62;
t98 = t58 * t105;
t33 = -t79 * t103 + t64 * t98;
t56 = pkin(13) + qJ(5);
t53 = sin(t56);
t54 = cos(t56);
t127 = t22 * t53 - t33 * t54;
t126 = t22 * t54 + t33 * t53;
t61 = sin(qJ(6));
t63 = cos(qJ(6));
t121 = m(7) * pkin(5) + t63 * mrSges(7,1) - t61 * mrSges(7,2) + mrSges(6,1);
t96 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t74 = t102 * t64 + t112 * t89;
t123 = t74 * t103 + t112 * t98;
t19 = t128 * t113 + t41 * t62;
t122 = -m(6) - m(7);
t120 = t74 * t105 - t112 * t97;
t85 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t116 = -t61 * mrSges(7,1) - t63 * mrSges(7,2) + t85;
t118 = m(5) - t122;
t57 = sin(pkin(13));
t59 = cos(pkin(13));
t117 = m(5) * pkin(3) + t59 * mrSges(5,1) - t57 * mrSges(5,2) + t121 * t54 - t96 * t53 + mrSges(4,1);
t111 = t33 * t57;
t110 = t58 * t64;
t100 = t58 * t112;
t107 = t64 * pkin(1) + qJ(2) * t100;
t93 = -pkin(1) * t112 + qJ(2) * t110;
t88 = t106 * t103;
t86 = t105 * t104;
t70 = -t41 * pkin(2) + t33 * pkin(9) + t93;
t42 = t104 * t64 - t112 * t87;
t69 = t42 * pkin(2) + t123 * pkin(9) + t107;
t23 = t120 * t113 + t42 * t62;
t24 = t42 * t113 - t120 * t62;
t52 = pkin(4) * t59 + pkin(3);
t60 = -pkin(10) - qJ(4);
t65 = t123 * t57;
t66 = pkin(4) * t65 - t23 * t60 + t24 * t52 + t69;
t40 = -t104 * t97 + t105 * t106;
t32 = t62 * t88 + (t102 * t113 + t62 * t86) * t58;
t31 = t102 * t58 * t62 + (-t58 * t86 - t88) * t113;
t14 = t32 * t54 + t40 * t53;
t8 = t123 * t53 + t24 * t54;
t7 = -t123 * t54 + t24 * t53;
t2 = t23 * t61 + t63 * t8;
t1 = t23 * t63 - t61 * t8;
t3 = [(-t64 * mrSges(2,1) + t112 * mrSges(2,2) - m(3) * t107 - t42 * mrSges(3,1) + t74 * mrSges(3,2) - mrSges(3,3) * t100 - m(4) * t69 - t24 * mrSges(4,1) - t123 * mrSges(4,3) - m(5) * (t24 * pkin(3) + t69) - (t24 * t59 + t65) * mrSges(5,1) - (t123 * t59 - t24 * t57) * mrSges(5,2) - m(6) * t66 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t66) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t96 * t7 + t85 * t23) * g(2) + (t112 * mrSges(2,1) + t64 * mrSges(2,2) - m(3) * t93 + t41 * mrSges(3,1) - t79 * mrSges(3,2) - mrSges(3,3) * t110 - m(4) * t70 - t22 * mrSges(4,1) - t33 * mrSges(4,3) - m(5) * (t22 * pkin(3) + t70) - (t22 * t59 + t111) * mrSges(5,1) - (-t22 * t57 + t33 * t59) * mrSges(5,2) + t96 * t127 - t121 * t126 - t116 * t19 + t122 * (pkin(4) * t111 + t19 * t60 + t22 * t52 + t70)) * g(1) (-t106 * g(3) + (-t112 * g(1) + t64 * g(2)) * t58) * (m(3) + m(4) + t118) (t122 * (-t31 * t52 - t32 * t60) + t116 * t32 + t117 * t31) * g(3) + (t122 * (-t19 * t52 + t22 * t60) - t116 * t22 + t117 * t19) * g(2) + (t122 * (-t23 * t52 - t24 * t60) + t116 * t24 + t117 * t23) * g(1), t118 * (-g(1) * t23 - g(2) * t19 - g(3) * t31) (t96 * t14 - t121 * (-t32 * t53 + t40 * t54)) * g(3) + (-t121 * t127 - t126 * t96) * g(2) + (t121 * t7 + t96 * t8) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t126 * t61 + t19 * t63) * mrSges(7,1) + (t126 * t63 - t19 * t61) * mrSges(7,2)) - g(3) * ((-t14 * t61 + t31 * t63) * mrSges(7,1) + (-t14 * t63 - t31 * t61) * mrSges(7,2))];
taug  = t3(:);
