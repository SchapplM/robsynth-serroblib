% Calculate joint inertia matrix for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2018-11-23 15:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:53:11
% EndTime: 2018-11-23 15:53:12
% DurationCPUTime: 0.44s
% Computational Cost: add. (448->160), mult. (756->196), div. (0->0), fcn. (507->6), ass. (0->60)
t41 = sin(qJ(6));
t43 = cos(qJ(6));
t56 = t41 ^ 2 + t43 ^ 2;
t12 = m(7) * t56;
t68 = m(6) + t12;
t76 = m(5) + t68;
t42 = sin(qJ(3));
t35 = t42 ^ 2;
t44 = cos(qJ(3));
t37 = t44 ^ 2;
t75 = t37 + t35;
t53 = t56 * mrSges(7,3);
t74 = t53 - mrSges(6,2);
t73 = 0.2e1 * t75;
t38 = sin(pkin(9));
t23 = t38 * pkin(1) + pkin(7);
t72 = -qJ(5) + t23;
t39 = cos(pkin(9));
t24 = -t39 * pkin(1) - pkin(2);
t25 = t42 * qJ(4);
t57 = t44 * pkin(3) + t25;
t7 = -t57 + t24;
t71 = -0.2e1 * t7;
t10 = t72 * t44;
t70 = t10 ^ 2;
t69 = -0.2e1 * t10;
t45 = -pkin(3) - pkin(4);
t67 = Ifges(7,4) * t41;
t66 = Ifges(7,4) * t43;
t65 = Ifges(7,5) * t41;
t64 = Ifges(7,5) * t43;
t63 = Ifges(7,6) * t43;
t62 = t10 * t42;
t61 = t44 * mrSges(7,3);
t60 = -mrSges(4,2) + mrSges(5,3);
t59 = t75 * t23 ^ 2;
t58 = Ifges(7,6) * t41 * t44 + Ifges(7,3) * t42;
t54 = -m(5) * pkin(3) - mrSges(5,1);
t51 = m(5) * t23 + mrSges(5,2) - mrSges(6,3);
t4 = t44 * pkin(4) - t7;
t3 = t42 * pkin(5) + t44 * pkin(8) + t4;
t9 = t72 * t42;
t1 = t43 * t3 - t41 * t9;
t2 = t41 * t3 + t43 * t9;
t50 = t41 * t1 - t43 * t2;
t49 = t41 * mrSges(7,1) + t43 * mrSges(7,2);
t13 = -t42 * mrSges(7,2) + t41 * t61;
t14 = t42 * mrSges(7,1) + t43 * t61;
t48 = -t43 * t13 + t41 * t14;
t46 = qJ(4) ^ 2;
t40 = qJ(4) + pkin(5);
t33 = -pkin(8) + t45;
t28 = t42 * mrSges(6,1);
t17 = -Ifges(7,1) * t41 - t66;
t16 = -Ifges(7,2) * t43 - t67;
t15 = t43 * mrSges(7,1) - t41 * mrSges(7,2);
t8 = t49 * t44;
t6 = Ifges(7,5) * t42 + (-Ifges(7,1) * t43 + t67) * t44;
t5 = Ifges(7,6) * t42 + (Ifges(7,2) * t41 - t66) * t44;
t11 = [0.2e1 * t1 * t14 + t8 * t69 + 0.2e1 * t2 * t13 + 0.2e1 * t4 * t28 + Ifges(2,3) + Ifges(3,3) + (-0.2e1 * t24 * mrSges(4,1) + mrSges(5,1) * t71 - 0.2e1 * t4 * mrSges(6,2) + mrSges(6,3) * t69 + t41 * t5 - t43 * t6 + (Ifges(6,1) + Ifges(5,3) + Ifges(4,2)) * t44) * t44 + (0.2e1 * t24 * mrSges(4,2) + mrSges(5,3) * t71 - 0.2e1 * t9 * mrSges(6,3) + (Ifges(6,2) + Ifges(5,1) + Ifges(4,1)) * t42 + ((2 * Ifges(4,4)) + (2 * Ifges(6,4)) - (2 * Ifges(5,5)) - t64) * t44 + t58) * t42 + m(7) * (t1 ^ 2 + t2 ^ 2 + t70) + m(6) * (t4 ^ 2 + t9 ^ 2 + t70) + m(5) * (t7 ^ 2 + t59) + m(4) * (t24 ^ 2 + t59) + (mrSges(5,2) + mrSges(4,3)) * t23 * t73 + (0.2e1 * t39 * mrSges(3,1) - 0.2e1 * t38 * mrSges(3,2) + m(3) * (t38 ^ 2 + t39 ^ 2) * pkin(1)) * pkin(1); -t42 * t8 + t48 * t44 + m(6) * (-t9 * t44 + t62) + m(7) * (t44 * t50 + t62); m(3) + m(7) * (t37 * t56 + t35) + (m(6) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t73; t9 * mrSges(6,2) - t40 * t8 + (t15 + mrSges(6,1)) * t10 + (-t5 / 0.2e1 - t2 * mrSges(7,3) + t33 * t13) * t43 + (-t6 / 0.2e1 + t1 * mrSges(7,3) - t33 * t14) * t41 + m(7) * (t40 * t10 - t33 * t50) + m(6) * (qJ(4) * t10 + t45 * t9) + (-t65 / 0.2e1 - t63 / 0.2e1 + Ifges(6,6) + Ifges(5,4) + Ifges(4,5) - pkin(3) * mrSges(5,2) - t45 * mrSges(6,3) + (-mrSges(4,1) + t54) * t23) * t42 + (Ifges(6,5) - Ifges(5,6) + Ifges(4,6) + t41 * t16 / 0.2e1 - t43 * t17 / 0.2e1 + t60 * t23 + t51 * qJ(4)) * t44; t28 + (t15 + t60) * t42 + (mrSges(4,1) + mrSges(5,1) + t74) * t44 + m(6) * (-t45 * t44 + t25) + m(7) * (-t33 * t44 * t56 + t40 * t42) + m(5) * t57; 0.2e1 * pkin(3) * mrSges(5,1) + 0.2e1 * t45 * mrSges(6,2) + 0.2e1 * t40 * t15 - t43 * t16 - t41 * t17 + Ifges(5,2) + Ifges(4,3) + Ifges(6,3) + m(7) * (t33 ^ 2 * t56 + t40 ^ 2) + m(5) * (pkin(3) ^ 2 + t46) + m(6) * (t45 ^ 2 + t46) + 0.2e1 * (mrSges(5,3) + mrSges(6,1)) * qJ(4) - 0.2e1 * t33 * t53; m(6) * t9 - m(7) * t50 + t42 * t51 - t48; -t76 * t44; m(6) * t45 + t12 * t33 + t54 - t74; t76; -t44 * mrSges(6,2) + t41 * t13 + t43 * t14 + t28 + m(7) * (t43 * t1 + t41 * t2) + m(6) * t4; 0; 0; 0; t68; t1 * mrSges(7,1) - t2 * mrSges(7,2) - t44 * t64 + t58; t8; -t33 * t49 - t63 - t65; -t49; t15; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1) t11(2) t11(4) t11(7) t11(11) t11(16); t11(2) t11(3) t11(5) t11(8) t11(12) t11(17); t11(4) t11(5) t11(6) t11(9) t11(13) t11(18); t11(7) t11(8) t11(9) t11(10) t11(14) t11(19); t11(11) t11(12) t11(13) t11(14) t11(15) t11(20); t11(16) t11(17) t11(18) t11(19) t11(20) t11(21);];
Mq  = res;
