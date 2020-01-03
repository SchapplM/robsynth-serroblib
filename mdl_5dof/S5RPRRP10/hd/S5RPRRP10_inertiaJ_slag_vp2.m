% Calculate joint inertia matrix for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:56
% EndTime: 2019-12-31 18:50:58
% DurationCPUTime: 0.68s
% Computational Cost: add. (627->159), mult. (1218->218), div. (0->0), fcn. (1206->6), ass. (0->68)
t94 = Ifges(5,5) + Ifges(6,5);
t93 = Ifges(5,6) + Ifges(6,6);
t92 = Ifges(5,3) + Ifges(6,3);
t91 = -2 * mrSges(6,3);
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t68 = t58 ^ 2 + t60 ^ 2;
t90 = 0.2e1 * t68;
t89 = t94 * t58 + t93 * t60;
t57 = cos(pkin(8));
t73 = pkin(6) + qJ(2);
t33 = t73 * t57;
t59 = sin(qJ(3));
t56 = sin(pkin(8));
t66 = t73 * t56;
t83 = cos(qJ(3));
t20 = t59 * t33 + t83 * t66;
t88 = t20 ^ 2;
t53 = t57 ^ 2;
t87 = 0.2e1 * t20;
t44 = -t57 * pkin(2) - pkin(1);
t86 = 0.2e1 * t44;
t85 = m(6) * pkin(4);
t82 = Ifges(5,4) * t58;
t81 = Ifges(5,4) * t60;
t80 = Ifges(6,4) * t58;
t79 = Ifges(6,4) * t60;
t32 = t83 * t56 + t59 * t57;
t78 = t32 * t58;
t77 = t32 * t60;
t76 = t58 * mrSges(5,2);
t75 = t58 * mrSges(6,3);
t72 = -qJ(5) - pkin(7);
t31 = t59 * t56 - t83 * t57;
t14 = t31 * pkin(3) - t32 * pkin(7) + t44;
t22 = t83 * t33 - t59 * t66;
t4 = t58 * t14 + t60 * t22;
t12 = mrSges(6,1) * t78 + mrSges(6,2) * t77;
t69 = t56 ^ 2 + t53;
t67 = qJ(5) * t32;
t65 = -t57 * mrSges(3,1) + t56 * mrSges(3,2);
t47 = t58 * mrSges(6,2);
t35 = -t60 * mrSges(6,1) + t47;
t3 = t60 * t14 - t58 * t22;
t64 = t92 * t31 + t94 * t77;
t63 = mrSges(5,1) * t58 + mrSges(5,2) * t60;
t45 = -t60 * pkin(4) - pkin(3);
t41 = Ifges(5,1) * t58 + t81;
t40 = Ifges(6,1) * t58 + t79;
t39 = Ifges(5,2) * t60 + t82;
t38 = Ifges(6,2) * t60 + t80;
t37 = t72 * t60;
t36 = -t60 * mrSges(5,1) + t76;
t34 = t72 * t58;
t27 = t32 * mrSges(4,2);
t18 = t31 * mrSges(5,1) - mrSges(5,3) * t77;
t17 = t31 * mrSges(6,1) - mrSges(6,3) * t77;
t16 = -t31 * mrSges(5,2) - mrSges(5,3) * t78;
t15 = -t31 * mrSges(6,2) - t32 * t75;
t13 = t63 * t32;
t9 = pkin(4) * t78 + t20;
t8 = Ifges(5,5) * t31 + (Ifges(5,1) * t60 - t82) * t32;
t7 = Ifges(6,5) * t31 + (Ifges(6,1) * t60 - t80) * t32;
t6 = Ifges(5,6) * t31 + (-Ifges(5,2) * t58 + t81) * t32;
t5 = Ifges(6,6) * t31 + (-Ifges(6,2) * t58 + t79) * t32;
t2 = -t58 * t67 + t4;
t1 = t31 * pkin(4) - t60 * t67 + t3;
t10 = [-0.2e1 * pkin(1) * t65 + Ifges(3,2) * t53 + t27 * t86 + 0.2e1 * t9 * t12 + 0.2e1 * t2 * t15 + 0.2e1 * t4 * t16 + 0.2e1 * t1 * t17 + 0.2e1 * t3 * t18 + t13 * t87 + Ifges(2,3) + (Ifges(3,1) * t56 + 0.2e1 * Ifges(3,4) * t57) * t56 + 0.2e1 * t69 * qJ(2) * mrSges(3,3) + (mrSges(4,1) * t86 - 0.2e1 * t22 * mrSges(4,3) + Ifges(4,2) * t31 + t64) * t31 + m(6) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2 + t88) + m(4) * (t22 ^ 2 + t44 ^ 2 + t88) + m(3) * (t69 * qJ(2) ^ 2 + pkin(1) ^ 2) + (mrSges(4,3) * t87 + Ifges(4,1) * t32 - 0.2e1 * Ifges(4,4) * t31 + (t7 + t8) * t60 + (-t31 * t93 - t5 - t6) * t58) * t32; -m(3) * pkin(1) + t31 * mrSges(4,1) + t27 + (t17 + t18) * t60 + (t15 + t16) * t58 + m(5) * (t60 * t3 + t58 * t4) + m(6) * (t60 * t1 + t58 * t2) + m(4) * t44 + t65; m(3) + m(4) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t90; t45 * t12 + t34 * t17 + t9 * t35 - t37 * t15 - pkin(3) * t13 - t22 * mrSges(4,2) + m(6) * (t34 * t1 - t37 * t2 + t45 * t9) + (-m(5) * pkin(3) - mrSges(4,1) + t36) * t20 + (t5 / 0.2e1 + t6 / 0.2e1 + t2 * mrSges(6,3) + t4 * mrSges(5,3) + (m(5) * t4 + t16) * pkin(7)) * t60 + (t7 / 0.2e1 + t8 / 0.2e1 - t1 * mrSges(6,3) - t3 * mrSges(5,3) + (-m(5) * t3 - t18) * pkin(7)) * t58 + (Ifges(4,5) + (t40 / 0.2e1 + t41 / 0.2e1) * t60 + (-t38 / 0.2e1 - t39 / 0.2e1) * t58) * t32 + (-Ifges(4,6) + t89 / 0.2e1) * t31; m(6) * (t34 * t60 - t37 * t58); -0.2e1 * pkin(3) * t36 + 0.2e1 * t45 * t35 + Ifges(4,3) + pkin(7) * mrSges(5,3) * t90 + m(6) * (t34 ^ 2 + t37 ^ 2 + t45 ^ 2) + m(5) * (t68 * pkin(7) ^ 2 + pkin(3) ^ 2) + (t37 * t91 + t38 + t39) * t60 + (t34 * t91 + t40 + t41) * t58; t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) - t93 * t78 + (m(6) * t1 + t17) * pkin(4) + t64; -t76 - t47 + (mrSges(5,1) + mrSges(6,1) + t85) * t60; t34 * mrSges(6,1) + t37 * mrSges(6,2) - t63 * pkin(7) + (m(6) * t34 - t75) * pkin(4) + t89; (0.2e1 * mrSges(6,1) + t85) * pkin(4) + t92; m(6) * t9 + t12; 0; m(6) * t45 + t35; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
