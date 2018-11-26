% Calculate joint inertia matrix for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:42
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:42:04
% EndTime: 2018-11-23 16:42:04
% DurationCPUTime: 0.97s
% Computational Cost: add. (1379->247), mult. (2631->346), div. (0->0), fcn. (2782->8), ass. (0->94)
t83 = sin(pkin(10));
t85 = cos(pkin(10));
t122 = t83 ^ 2 + t85 ^ 2;
t121 = 0.2e1 * t122;
t106 = -qJ(3) - pkin(7);
t90 = cos(qJ(2));
t71 = t106 * t90;
t84 = sin(pkin(9));
t86 = cos(pkin(9));
t88 = sin(qJ(2));
t99 = t106 * t88;
t40 = -t71 * t84 - t86 * t99;
t120 = t40 ^ 2;
t119 = 0.2e1 * t40;
t77 = -pkin(2) * t90 - pkin(1);
t118 = 0.2e1 * t77;
t117 = pkin(4) + pkin(5);
t116 = pkin(2) * t84;
t115 = pkin(2) * t86;
t74 = qJ(4) + t116;
t114 = -pkin(8) + t74;
t113 = Ifges(5,4) * t83;
t112 = Ifges(5,4) * t85;
t111 = Ifges(6,5) * t83;
t110 = Ifges(6,5) * t85;
t60 = t84 * t90 + t86 * t88;
t109 = t60 * t83;
t108 = t60 * t85;
t107 = t83 * mrSges(6,3);
t58 = t84 * t88 - t86 * t90;
t31 = -mrSges(5,2) * t58 - mrSges(5,3) * t109;
t34 = -mrSges(6,2) * t109 + mrSges(6,3) * t58;
t105 = t31 + t34;
t32 = mrSges(5,1) * t58 - mrSges(5,3) * t108;
t33 = -t58 * mrSges(6,1) + mrSges(6,2) * t108;
t104 = t32 - t33;
t30 = pkin(3) * t58 - qJ(4) * t60 + t77;
t42 = -t86 * t71 + t84 * t99;
t13 = t83 * t30 + t85 * t42;
t29 = mrSges(5,1) * t109 + mrSges(5,2) * t108;
t103 = t122 * t74 ^ 2;
t101 = t88 ^ 2 + t90 ^ 2;
t7 = t58 * qJ(5) + t13;
t87 = sin(qJ(6));
t89 = cos(qJ(6));
t61 = t83 * t89 - t85 * t87;
t23 = t61 * t60;
t94 = t83 * t87 + t85 * t89;
t24 = t94 * t60;
t100 = -Ifges(7,5) * t24 - Ifges(7,6) * t23 + Ifges(7,3) * t58;
t76 = -pkin(3) - t115;
t35 = mrSges(7,1) * t94 + t61 * mrSges(7,2);
t38 = t83 * t42;
t12 = t30 * t85 - t38;
t28 = mrSges(6,1) * t109 - mrSges(6,3) * t108;
t97 = -qJ(5) * t108 + t40;
t8 = -pkin(4) * t58 - t12;
t96 = t7 * t85 + t8 * t83;
t9 = -t23 * mrSges(7,1) + t24 * mrSges(7,2);
t95 = -t12 * t83 + t13 * t85;
t93 = qJ(5) * t83 - t76;
t78 = t83 * mrSges(5,2);
t68 = Ifges(5,1) * t83 + t112;
t67 = Ifges(6,1) * t83 - t110;
t66 = Ifges(5,2) * t85 + t113;
t65 = -Ifges(6,3) * t85 + t111;
t64 = -t85 * mrSges(5,1) + t78;
t63 = -t85 * mrSges(6,1) - t107;
t57 = Ifges(7,5) * t61;
t56 = Ifges(7,6) * t94;
t54 = t60 * mrSges(4,2);
t51 = t114 * t85;
t50 = t114 * t83;
t49 = -pkin(4) * t85 - t93;
t43 = t117 * t85 + t93;
t37 = Ifges(7,1) * t61 - Ifges(7,4) * t94;
t36 = Ifges(7,4) * t61 - Ifges(7,2) * t94;
t27 = t50 * t87 + t51 * t89;
t26 = t50 * t89 - t51 * t87;
t20 = Ifges(5,5) * t58 + (Ifges(5,1) * t85 - t113) * t60;
t19 = Ifges(6,4) * t58 + (Ifges(6,1) * t85 + t111) * t60;
t18 = Ifges(5,6) * t58 + (-Ifges(5,2) * t83 + t112) * t60;
t17 = Ifges(6,6) * t58 + (Ifges(6,3) * t83 + t110) * t60;
t16 = -mrSges(7,1) * t58 - mrSges(7,3) * t24;
t15 = mrSges(7,2) * t58 + mrSges(7,3) * t23;
t14 = pkin(4) * t109 + t97;
t10 = t109 * t117 + t97;
t6 = Ifges(7,1) * t24 + Ifges(7,4) * t23 - Ifges(7,5) * t58;
t5 = Ifges(7,4) * t24 + Ifges(7,2) * t23 - Ifges(7,6) * t58;
t4 = pkin(8) * t109 + t7;
t3 = t38 + (-pkin(8) * t60 - t30) * t85 - t117 * t58;
t2 = t3 * t87 + t4 * t89;
t1 = t3 * t89 - t4 * t87;
t11 = [t90 * (Ifges(3,4) * t88 + Ifges(3,2) * t90) - 0.2e1 * pkin(1) * (-t90 * mrSges(3,1) + t88 * mrSges(3,2)) + t88 * (Ifges(3,1) * t88 + Ifges(3,4) * t90) + t54 * t118 + t29 * t119 + t23 * t5 + t24 * t6 + 0.2e1 * t14 * t28 + 0.2e1 * t13 * t31 + 0.2e1 * t12 * t32 + 0.2e1 * t8 * t33 + 0.2e1 * t7 * t34 + 0.2e1 * t2 * t15 + 0.2e1 * t1 * t16 - 0.2e1 * t10 * t9 + Ifges(2,3) + 0.2e1 * t101 * pkin(7) * mrSges(3,3) + (mrSges(4,1) * t118 - 0.2e1 * t42 * mrSges(4,3) + (Ifges(5,3) + Ifges(6,2) + Ifges(4,2)) * t58 + t100) * t58 + m(3) * (t101 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t42 ^ 2 + t77 ^ 2 + t120) + m(5) * (t12 ^ 2 + t13 ^ 2 + t120) + m(6) * (t14 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(7) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + (mrSges(4,3) * t119 + Ifges(4,1) * t60 + (t19 + t20) * t85 + (t17 - t18) * t83 + (-(2 * Ifges(4,4)) + (Ifges(6,4) + Ifges(5,5)) * t85 + (-Ifges(5,6) + Ifges(6,6)) * t83) * t58) * t60; m(7) * (t1 * t26 - t10 * t43 + t2 * t27) + (-mrSges(4,3) * t115 + Ifges(4,5) + (t67 / 0.2e1 + t68 / 0.2e1) * t85 + (t65 / 0.2e1 - t66 / 0.2e1) * t83) * t60 + (-t104 * t83 + t105 * t85) * t74 + t96 * mrSges(6,2) + m(6) * (t14 * t49 + t74 * t96) + (-t88 * mrSges(3,1) - t90 * mrSges(3,2)) * pkin(7) + (t19 / 0.2e1 + t20 / 0.2e1) * t83 + (-t17 / 0.2e1 + t18 / 0.2e1) * t85 + t95 * mrSges(5,3) + m(5) * (t40 * t76 + t74 * t95) + (-t1 * t61 - t2 * t94) * mrSges(7,3) - t94 * t5 / 0.2e1 + (-t57 / 0.2e1 + t56 / 0.2e1 - Ifges(4,6) - mrSges(4,3) * t116 + (Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t85 + (Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t83) * t58 + Ifges(3,6) * t90 + Ifges(3,5) * t88 + t76 * t29 + t61 * t6 / 0.2e1 + t14 * t63 - t10 * t35 + t23 * t36 / 0.2e1 + (t64 - mrSges(4,1)) * t40 + t24 * t37 / 0.2e1 - t42 * mrSges(4,2) + t43 * t9 + t49 * t28 + t26 * t16 + t27 * t15 + m(4) * (-t40 * t86 + t42 * t84) * pkin(2); 0.2e1 * t43 * t35 - t94 * t36 + t61 * t37 + 0.2e1 * t49 * t63 + 0.2e1 * t76 * t64 + Ifges(3,3) + Ifges(4,3) + (t66 - t65) * t85 + (t68 + t67) * t83 + m(7) * (t26 ^ 2 + t27 ^ 2 + t43 ^ 2) + m(5) * (t76 ^ 2 + t103) + m(6) * (t49 ^ 2 + t103) + m(4) * (t84 ^ 2 + t86 ^ 2) * pkin(2) ^ 2 + (mrSges(6,2) + mrSges(5,3)) * t74 * t121 + 0.2e1 * (mrSges(4,1) * t86 - mrSges(4,2) * t84) * pkin(2) + 0.2e1 * (-t26 * t61 - t27 * t94) * mrSges(7,3); t58 * mrSges(4,1) + t61 * t15 - t94 * t16 + t54 + t104 * t85 + t105 * t83 + m(7) * (-t1 * t94 + t2 * t61) + m(6) * (t7 * t83 - t8 * t85) + m(5) * (t12 * t85 + t13 * t83) + m(4) * t77; m(7) * (-t26 * t94 + t27 * t61); m(4) + m(7) * (t61 ^ 2 + t94 ^ 2) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t121; m(5) * t40 + m(6) * t14 + m(7) * t10 + t28 + t29 - t9; -t107 + t78 + (-mrSges(6,1) - mrSges(5,1)) * t85 + m(5) * t76 + m(6) * t49 - m(7) * t43 - t35; 0; m(5) + m(6) + m(7); t87 * t15 + t89 * t16 + m(7) * (t1 * t89 + t2 * t87) + m(6) * t8 + t33; m(7) * (t26 * t89 + t27 * t87) + (m(6) * t74 + mrSges(6,2)) * t83 + (-t61 * t89 - t87 * t94) * mrSges(7,3); -m(6) * t85 + m(7) * (t61 * t87 - t89 * t94); 0; m(6) + m(7) * (t87 ^ 2 + t89 ^ 2); mrSges(7,1) * t1 - mrSges(7,2) * t2 - t100; mrSges(7,1) * t26 - mrSges(7,2) * t27 - t56 + t57; -t35; 0; mrSges(7,1) * t89 - t87 * mrSges(7,2); Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1) t11(2) t11(4) t11(7) t11(11) t11(16); t11(2) t11(3) t11(5) t11(8) t11(12) t11(17); t11(4) t11(5) t11(6) t11(9) t11(13) t11(18); t11(7) t11(8) t11(9) t11(10) t11(14) t11(19); t11(11) t11(12) t11(13) t11(14) t11(15) t11(20); t11(16) t11(17) t11(18) t11(19) t11(20) t11(21);];
Mq  = res;
