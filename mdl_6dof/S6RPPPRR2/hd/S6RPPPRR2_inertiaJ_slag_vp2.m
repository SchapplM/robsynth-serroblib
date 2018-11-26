% Calculate joint inertia matrix for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2018-11-23 15:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:37:41
% EndTime: 2018-11-23 15:37:42
% DurationCPUTime: 0.38s
% Computational Cost: add. (663->135), mult. (1154->190), div. (0->0), fcn. (1098->8), ass. (0->61)
t43 = sin(pkin(10));
t45 = cos(pkin(10));
t48 = sin(qJ(5));
t66 = t48 * t45;
t72 = cos(qJ(5));
t20 = t43 * t72 + t66;
t59 = t72 * t45;
t19 = t48 * t43 - t59;
t47 = sin(qJ(6));
t69 = t19 * t47;
t10 = -t20 * mrSges(7,2) + mrSges(7,3) * t69;
t49 = cos(qJ(6));
t68 = t19 * t49;
t11 = t20 * mrSges(7,1) + mrSges(7,3) * t68;
t52 = -t49 * t10 + t47 * t11;
t25 = -t49 * mrSges(7,1) + t47 * mrSges(7,2);
t84 = -m(7) * pkin(5) - mrSges(6,1) + t25;
t17 = t20 ^ 2;
t18 = t19 ^ 2;
t40 = t45 ^ 2;
t61 = t43 ^ 2 + t40;
t23 = m(5) * t61;
t83 = m(4) + t23 + m(6) * (t17 + t18);
t46 = cos(pkin(9));
t33 = -t46 * pkin(1) - pkin(2);
t29 = -qJ(4) + t33;
t73 = -pkin(7) + t29;
t14 = t73 * t43;
t6 = t48 * t14 - t59 * t73;
t82 = t6 ^ 2;
t44 = sin(pkin(9));
t30 = t44 * pkin(1) + qJ(3);
t81 = t30 ^ 2;
t80 = -2 * mrSges(6,3);
t78 = t47 / 0.2e1;
t77 = pkin(5) * t20;
t76 = t19 * t6;
t75 = t6 * t20;
t71 = Ifges(7,4) * t47;
t70 = Ifges(7,4) * t49;
t64 = -Ifges(7,5) * t68 + Ifges(7,3) * t20;
t63 = t43 * mrSges(5,1) + t45 * mrSges(5,2);
t62 = Ifges(7,5) * t47 + Ifges(7,6) * t49;
t60 = t47 ^ 2 + t49 ^ 2;
t58 = t61 * mrSges(5,3);
t57 = t60 * t19;
t15 = t19 * mrSges(6,2);
t55 = -t20 * mrSges(6,1) + t15;
t24 = t43 * pkin(4) + t30;
t5 = t19 * pkin(8) + t24 + t77;
t8 = t14 * t72 + t66 * t73;
t1 = -t47 * t8 + t49 * t5;
t2 = t47 * t5 + t49 * t8;
t54 = t1 * t47 - t2 * t49;
t53 = -mrSges(7,1) * t47 - mrSges(7,2) * t49;
t27 = Ifges(7,1) * t47 + t70;
t26 = Ifges(7,2) * t49 + t71;
t9 = t53 * t19;
t4 = Ifges(7,5) * t20 + (-Ifges(7,1) * t49 + t71) * t19;
t3 = Ifges(7,6) * t20 + (Ifges(7,2) * t47 - t70) * t19;
t7 = [Ifges(5,1) * t40 + 0.2e1 * t33 * mrSges(4,2) - 0.2e1 * t24 * t15 + 0.2e1 * t6 * t9 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (-0.2e1 * Ifges(5,4) * t45 + Ifges(5,2) * t43) * t43 + (0.2e1 * t24 * mrSges(6,1) + Ifges(6,2) * t20 + t8 * t80 + t64) * t20 + (t6 * t80 + Ifges(6,1) * t19 + t47 * t3 - t49 * t4 + (Ifges(7,6) * t47 + (2 * Ifges(6,4))) * t20) * t19 + m(7) * (t1 ^ 2 + t2 ^ 2 + t82) + m(6) * (t24 ^ 2 + t8 ^ 2 + t82) + m(4) * (t33 ^ 2 + t81) + m(5) * (t29 ^ 2 * t61 + t81) + m(3) * (t44 ^ 2 + t46 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(4,3) + t63) * t30 + 0.2e1 * (t46 * mrSges(3,1) - t44 * mrSges(3,2)) * pkin(1) - 0.2e1 * t29 * t58; t20 * t9 + t52 * t19 + m(7) * (t19 * t54 + t75) + m(6) * (-t8 * t19 + t75); m(3) + m(7) * (t18 * t60 + t17) + t83; -t18 * mrSges(6,3) + t19 * t9 + mrSges(4,2) - t58 + (-mrSges(6,3) * t20 - t52) * t20 + m(7) * (-t20 * t54 + t76) + m(6) * (t20 * t8 + t76) + m(4) * t33 + t29 * t23; m(7) * (0.1e1 - t60) * t20 * t19; m(7) * (t17 * t60 + t18) + t83; t47 * t10 + t49 * t11 + m(7) * (t49 * t1 + t47 * t2) + m(6) * t24 + m(5) * t30 - t55 + t63; 0; 0; m(7) * t60 + m(5) + m(6); t4 * t78 + t49 * t3 / 0.2e1 - pkin(5) * t9 - t8 * mrSges(6,2) - t54 * mrSges(7,3) + (-t49 * t27 / 0.2e1 + t26 * t78 - Ifges(6,5)) * t19 + (t62 / 0.2e1 - Ifges(6,6)) * t20 + t84 * t6 + (-m(7) * t54 - t52) * pkin(8); m(7) * (-pkin(8) * t57 - t77) + t20 * t25 - mrSges(7,3) * t57 + t55; t84 * t19 + (-mrSges(6,2) + (m(7) * pkin(8) + mrSges(7,3)) * t60) * t20; 0; Ifges(6,3) + t47 * t27 + t49 * t26 - 0.2e1 * pkin(5) * t25 + m(7) * (pkin(8) ^ 2 * t60 + pkin(5) ^ 2) + 0.2e1 * t60 * pkin(8) * mrSges(7,3); t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,6) * t69 + t64; -t9; t53 * t20; -t25; pkin(8) * t53 + t62; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
