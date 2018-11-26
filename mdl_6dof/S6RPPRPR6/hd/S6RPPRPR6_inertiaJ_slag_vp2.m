% Calculate joint inertia matrix for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2018-11-23 15:42
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:42:21
% EndTime: 2018-11-23 15:42:22
% DurationCPUTime: 0.40s
% Computational Cost: add. (355->140), mult. (586->169), div. (0->0), fcn. (358->4), ass. (0->58)
t35 = sin(qJ(4));
t27 = t35 ^ 2;
t37 = cos(qJ(4));
t29 = t37 ^ 2;
t49 = t29 + t27;
t74 = (mrSges(6,1) + mrSges(5,3)) * t49;
t34 = sin(qJ(6));
t36 = cos(qJ(6));
t50 = t34 ^ 2 + t36 ^ 2;
t8 = m(7) * t50;
t63 = m(6) + t8;
t73 = t50 * mrSges(7,3);
t72 = mrSges(5,1) - mrSges(6,2);
t70 = -mrSges(5,2) + mrSges(6,3);
t22 = qJ(5) * t35;
t69 = m(6) * (pkin(4) * t37 + t22);
t32 = pkin(1) + qJ(3);
t68 = t32 ^ 2;
t45 = -t37 * qJ(5) + t32;
t11 = t35 * pkin(4) + t45;
t67 = -0.2e1 * t11;
t66 = 0.2e1 * t32;
t65 = -t34 / 0.2e1;
t64 = t36 / 0.2e1;
t62 = pkin(4) + pkin(8);
t31 = qJ(2) - pkin(7);
t60 = pkin(5) - t31;
t59 = Ifges(7,4) * t34;
t58 = Ifges(7,4) * t36;
t55 = t35 * t36;
t10 = -t37 * mrSges(7,2) + mrSges(7,3) * t55;
t57 = t34 * t10;
t56 = t34 * t35;
t54 = t36 * t62;
t14 = t34 * mrSges(7,1) + t36 * mrSges(7,2);
t52 = t14 + mrSges(6,3);
t51 = t49 * t31 ^ 2;
t48 = Ifges(7,5) * t56 + Ifges(7,6) * t55 + Ifges(7,3) * t37;
t47 = t50 * t62;
t13 = t60 * t37;
t5 = t62 * t35 + t45;
t1 = t36 * t13 - t34 * t5;
t2 = t34 * t13 + t36 * t5;
t44 = t36 * t1 + t34 * t2;
t9 = t37 * mrSges(7,1) - mrSges(7,3) * t56;
t43 = -t36 * t9 - t57;
t42 = t36 * mrSges(7,1) - t34 * mrSges(7,2);
t41 = 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t49;
t40 = (qJ(2) ^ 2);
t39 = qJ(5) ^ 2;
t24 = Ifges(7,5) * t36;
t16 = Ifges(7,1) * t36 - t59;
t15 = -Ifges(7,2) * t34 + t58;
t12 = t60 * t35;
t6 = t42 * t35;
t4 = Ifges(7,5) * t37 + (Ifges(7,1) * t34 + t58) * t35;
t3 = Ifges(7,6) * t37 + (Ifges(7,2) * t36 + t59) * t35;
t7 = [-(2 * pkin(1) * mrSges(3,2)) + mrSges(4,3) * t66 + 0.2e1 * t1 * t9 + 0.2e1 * t2 * t10 + 0.2e1 * t12 * t6 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3) + (mrSges(5,2) * t66 + mrSges(6,3) * t67 + (Ifges(6,2) + Ifges(5,1)) * t37 + t48) * t37 + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) + m(6) * (t11 ^ 2 + t51) + m(5) * (t51 + t68) + (m(3) * (pkin(1) ^ 2 + t40)) + m(4) * (t40 + t68) + (2 * (mrSges(3,3) + mrSges(4,2)) * qJ(2)) - 0.2e1 * t31 * t74 + (mrSges(5,1) * t66 + mrSges(6,2) * t67 + t36 * t3 + t34 * t4 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t37 + (Ifges(6,3) + Ifges(5,2)) * t35) * t35; -(m(3) * pkin(1)) - t36 * t10 + t34 * t9 + mrSges(3,2) - mrSges(4,3) + t70 * t37 - t72 * t35 + m(7) * (t34 * t1 - t36 * t2) - m(6) * t11 + (-m(5) / 0.2e1 - m(4) / 0.2e1) * t66; m(3) + m(4) + m(5) + t63; m(4) * qJ(2) - t35 * t6 + mrSges(4,2) + t43 * t37 + m(7) * (-t35 * t12 - t44 * t37) + t31 * t41 - t74; 0; m(4) + m(7) * (t50 * t29 + t27) + t41; t4 * t64 + t3 * t65 - t12 * t14 - qJ(5) * t6 + m(7) * (-qJ(5) * t12 - t44 * t62) - t62 * t57 - t9 * t54 - t44 * mrSges(7,3) + (Ifges(7,6) * t65 + t24 / 0.2e1 + Ifges(5,5) - Ifges(6,4) - pkin(4) * mrSges(6,1)) * t37 + (-Ifges(5,6) + Ifges(6,5) + t34 * t16 / 0.2e1 + t15 * t64 - qJ(5) * mrSges(6,1)) * t35 + (t70 * t35 + t72 * t37 + t69) * t31; 0; (t72 + t73) * t37 + t69 + m(7) * (t37 * t47 + t22) + (-mrSges(5,2) + t52) * t35; -0.2e1 * pkin(4) * mrSges(6,2) - t34 * t15 + t36 * t16 + Ifges(6,1) + Ifges(5,3) + m(7) * (t50 * t62 ^ 2 + t39) + m(6) * (pkin(4) ^ 2 + t39) + 0.2e1 * mrSges(7,3) * t47 + 0.2e1 * t52 * qJ(5); m(7) * t44 + (-m(6) * t31 + mrSges(6,1)) * t37 - t43; 0; -t63 * t37; -m(6) * pkin(4) - t62 * t8 + mrSges(6,2) - t73; t63; t1 * mrSges(7,1) - t2 * mrSges(7,2) + t48; t14; -t42 * t37; -mrSges(7,1) * t54 + t24 + (mrSges(7,2) * t62 - Ifges(7,6)) * t34; t42; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
