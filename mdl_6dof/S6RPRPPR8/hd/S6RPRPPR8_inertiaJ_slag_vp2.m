% Calculate joint inertia matrix for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2018-11-23 15:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:56:02
% EndTime: 2018-11-23 15:56:03
% DurationCPUTime: 0.45s
% Computational Cost: add. (441->163), mult. (708->195), div. (0->0), fcn. (443->4), ass. (0->62)
t36 = sin(qJ(6));
t38 = cos(qJ(6));
t56 = t36 ^ 2 + t38 ^ 2;
t9 = m(7) * t56;
t70 = m(6) + t9;
t77 = m(5) + t70;
t37 = sin(qJ(3));
t31 = t37 ^ 2;
t39 = cos(qJ(3));
t33 = t39 ^ 2;
t55 = t33 + t31;
t49 = t56 * mrSges(7,3);
t76 = t49 - mrSges(6,2);
t41 = -pkin(1) - pkin(7);
t75 = t41 * t55;
t53 = -qJ(5) - t41;
t12 = t53 * t37;
t74 = t12 ^ 2;
t73 = 2 * mrSges(5,1);
t72 = -0.2e1 * t12;
t71 = 2 * qJ(2);
t40 = -pkin(3) - pkin(4);
t69 = Ifges(7,4) * t36;
t68 = Ifges(7,4) * t38;
t67 = Ifges(7,5) * t36;
t66 = Ifges(7,6) * t38;
t65 = t36 * t37;
t64 = t37 * t12;
t63 = t37 * t38;
t16 = t38 * mrSges(7,1) - t36 * mrSges(7,2);
t62 = mrSges(6,1) + t16;
t61 = mrSges(5,3) - mrSges(4,2);
t60 = mrSges(6,3) - mrSges(5,2);
t59 = Ifges(7,5) * t63 + Ifges(7,3) * t39;
t58 = t55 * t41 ^ 2;
t57 = t39 * mrSges(6,1) + t37 * mrSges(6,2);
t54 = t39 * qJ(4) - qJ(2);
t29 = -pkin(8) + t40;
t52 = m(5) / 0.2e1 + m(4) / 0.2e1;
t50 = -m(5) * pkin(3) - mrSges(5,1);
t47 = m(5) * t41 + t60;
t14 = t53 * t39;
t3 = t39 * pkin(5) + t29 * t37 + t54;
t1 = -t36 * t14 + t38 * t3;
t2 = t38 * t14 + t36 * t3;
t46 = t36 * t1 - t38 * t2;
t45 = t36 * mrSges(7,1) + t38 * mrSges(7,2);
t10 = -t39 * mrSges(7,2) - mrSges(7,3) * t65;
t11 = t39 * mrSges(7,1) - mrSges(7,3) * t63;
t44 = -t38 * t10 + t36 * t11;
t43 = qJ(2) ^ 2;
t42 = qJ(4) ^ 2;
t35 = qJ(4) + pkin(5);
t24 = qJ(4) * t37;
t18 = -Ifges(7,1) * t36 - t68;
t17 = -Ifges(7,2) * t38 - t69;
t15 = t37 * pkin(3) - t54;
t7 = t45 * t37;
t6 = t37 * t40 + t54;
t5 = Ifges(7,5) * t39 + (Ifges(7,1) * t38 - t69) * t37;
t4 = Ifges(7,6) * t39 + (-Ifges(7,2) * t36 + t68) * t37;
t8 = [0.2e1 * t6 * t57 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + t7 * t72 + Ifges(3,1) + Ifges(2,3) + (mrSges(3,3) * t71) - (2 * pkin(1) * mrSges(3,2)) + m(7) * (t1 ^ 2 + t2 ^ 2 + t74) + m(6) * (t14 ^ 2 + t6 ^ 2 + t74) + m(5) * (t15 ^ 2 + t58) + (m(3) * (pkin(1) ^ 2 + t43)) + m(4) * (t43 + t58) + ((mrSges(4,2) * t71) - 0.2e1 * t15 * mrSges(5,3) - 0.2e1 * t14 * mrSges(6,3) + (Ifges(6,2) + Ifges(5,1) + Ifges(4,1)) * t39 + t59) * t39 + ((mrSges(4,1) * t71) + t15 * t73 + mrSges(6,3) * t72 - t36 * t4 + t38 * t5 + (Ifges(6,1) + Ifges(5,3) + Ifges(4,2)) * t37 + (-Ifges(7,6) * t36 - (2 * Ifges(4,4)) - (2 * Ifges(6,4)) + (2 * Ifges(5,5))) * t39) * t37 - 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * t75; -(m(3) * pkin(1)) + t37 * t7 + mrSges(3,2) + t44 * t39 + m(7) * (t39 * t46 - t64) + m(6) * (-t39 * t14 - t64) + 0.2e1 * t52 * t75 + t55 * (-mrSges(4,3) + t60); m(3) + m(7) * (t33 * t56 + t31) + 0.2e1 * (m(6) / 0.2e1 + t52) * t55; t14 * mrSges(6,2) + t35 * t7 - t62 * t12 + (-t4 / 0.2e1 - t2 * mrSges(7,3) + t29 * t10) * t38 + (-t5 / 0.2e1 + t1 * mrSges(7,3) - t29 * t11) * t36 + m(7) * (-t35 * t12 - t46 * t29) + m(6) * (-qJ(4) * t12 + t40 * t14) + (-t67 / 0.2e1 - t66 / 0.2e1 + Ifges(6,6) + Ifges(5,4) + Ifges(4,5) - pkin(3) * mrSges(5,2) - t40 * mrSges(6,3) + (mrSges(4,1) - t50) * t41) * t39 + (-Ifges(6,5) + Ifges(5,6) - Ifges(4,6) - t36 * t17 / 0.2e1 + t38 * t18 / 0.2e1 + t61 * t41 + t47 * qJ(4)) * t37; (t61 + t62) * t37 + (mrSges(4,1) + mrSges(5,1) + t76) * t39 + m(6) * (-t40 * t39 + t24) + m(7) * (-t29 * t39 * t56 + t35 * t37) + m(5) * (pkin(3) * t39 + t24); pkin(3) * t73 + 0.2e1 * t40 * mrSges(6,2) + 0.2e1 * t35 * t16 - t38 * t17 - t36 * t18 + Ifges(5,2) + Ifges(4,3) + Ifges(6,3) + m(7) * (t29 ^ 2 * t56 + t35 ^ 2) + m(5) * (pkin(3) ^ 2 + t42) + m(6) * (t40 ^ 2 + t42) + 0.2e1 * (mrSges(5,3) + mrSges(6,1)) * qJ(4) - 0.2e1 * t29 * t49; m(6) * t14 - m(7) * t46 - t39 * t47 - t44; -t77 * t39; m(6) * t40 + t29 * t9 + t50 - t76; t77; t36 * t10 + t38 * t11 + m(7) * (t38 * t1 + t36 * t2) + m(6) * t6 + t57; 0; 0; 0; t70; t1 * mrSges(7,1) - t2 * mrSges(7,2) - Ifges(7,6) * t65 + t59; t45 * t39; -t29 * t45 - t66 - t67; -t45; t16; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t8(1) t8(2) t8(4) t8(7) t8(11) t8(16); t8(2) t8(3) t8(5) t8(8) t8(12) t8(17); t8(4) t8(5) t8(6) t8(9) t8(13) t8(18); t8(7) t8(8) t8(9) t8(10) t8(14) t8(19); t8(11) t8(12) t8(13) t8(14) t8(15) t8(20); t8(16) t8(17) t8(18) t8(19) t8(20) t8(21);];
Mq  = res;
