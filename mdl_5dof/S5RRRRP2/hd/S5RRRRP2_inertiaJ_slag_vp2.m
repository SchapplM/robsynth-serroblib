% Calculate joint inertia matrix for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m [6x1]
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:43
% EndTime: 2022-01-20 11:48:43
% DurationCPUTime: 0.43s
% Computational Cost: add. (657->148), mult. (1201->185), div. (0->0), fcn. (1082->6), ass. (0->60)
t60 = cos(qJ(3));
t96 = t60 ^ 2;
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t95 = (t59 * mrSges(5,1) + (-mrSges(5,2) - mrSges(6,2)) * t56) * pkin(3);
t58 = sin(qJ(2));
t46 = t58 * pkin(1) + pkin(7);
t57 = sin(qJ(3));
t32 = (-pkin(8) - t46) * t57;
t53 = t60 * pkin(8);
t33 = t60 * t46 + t53;
t12 = t59 * t32 - t56 * t33;
t13 = t56 * t32 + t59 * t33;
t37 = t56 * t60 + t59 * t57;
t76 = t37 * qJ(5);
t2 = t12 - t76;
t36 = -t56 * t57 + t59 * t60;
t26 = t36 * qJ(5);
t3 = t26 + t13;
t94 = t12 * mrSges(5,1) + t2 * mrSges(6,1) - t13 * mrSges(5,2) - t3 * mrSges(6,2);
t43 = (-pkin(8) - pkin(7)) * t57;
t44 = t60 * pkin(7) + t53;
t18 = t56 * t43 + t59 * t44;
t10 = t26 + t18;
t17 = t59 * t43 - t56 * t44;
t9 = t17 - t76;
t93 = t17 * mrSges(5,1) + t9 * mrSges(6,1) - t18 * mrSges(5,2) - t10 * mrSges(6,2);
t92 = 0.2e1 * mrSges(6,1);
t14 = -t36 * mrSges(6,1) + t37 * mrSges(6,2);
t91 = 0.2e1 * t14;
t15 = -t36 * mrSges(5,1) + t37 * mrSges(5,2);
t90 = 0.2e1 * t15;
t89 = m(5) * pkin(3);
t88 = m(6) * pkin(4);
t87 = pkin(3) * t56;
t85 = t59 * pkin(3);
t61 = cos(qJ(2));
t84 = t61 * pkin(1);
t80 = t37 * mrSges(6,3);
t78 = Ifges(5,3) + Ifges(6,3);
t77 = t57 ^ 2 + t96;
t75 = 2 * mrSges(5,3);
t74 = 0.2e1 * mrSges(6,3);
t49 = -t60 * pkin(3) - pkin(2);
t71 = t77 * t46;
t70 = (Ifges(5,5) + Ifges(6,5)) * t37 + (Ifges(5,6) + Ifges(6,6)) * t36;
t69 = -t57 * mrSges(4,1) - t60 * mrSges(4,2);
t68 = 0.2e1 * t77 * mrSges(4,3);
t20 = -t36 * pkin(4) + t49;
t67 = Ifges(4,2) * t96 + Ifges(3,3) + (Ifges(4,1) * t57 + 0.2e1 * Ifges(4,4) * t60) * t57 + (Ifges(6,1) + Ifges(5,1)) * t37 ^ 2 + ((Ifges(6,2) + Ifges(5,2)) * t36 + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t37) * t36;
t66 = (t61 * mrSges(3,1) - t58 * mrSges(3,2)) * pkin(1);
t47 = pkin(4) + t85;
t65 = Ifges(4,6) * t60 + Ifges(4,5) * t57 + (-mrSges(5,3) * t85 - t47 * mrSges(6,3)) * t37 + t70 + (mrSges(6,3) + mrSges(5,3)) * t36 * t87;
t63 = pkin(3) ^ 2;
t52 = t56 ^ 2 * t63;
t48 = -pkin(2) - t84;
t42 = -t60 * mrSges(4,1) + t57 * mrSges(4,2);
t41 = t49 - t84;
t19 = t20 - t84;
t1 = [t67 + m(4) * (t77 * t46 ^ 2 + t48 ^ 2) + t46 * t68 + (-t12 * t37 + t13 * t36) * t75 + (-t2 * t37 + t3 * t36) * t74 + m(3) * (t58 ^ 2 + t61 ^ 2) * pkin(1) ^ 2 + m(6) * (t19 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2 + t41 ^ 2) + 0.2e1 * t48 * t42 + t41 * t90 + t19 * t91 + Ifges(2,3) + 0.2e1 * t66; t67 + t66 + m(4) * (-pkin(2) * t48 + pkin(7) * t71) + (t77 * pkin(7) + t71) * mrSges(4,3) + ((-t12 - t17) * t37 + (t13 + t18) * t36) * mrSges(5,3) + ((-t2 - t9) * t37 + (t10 + t3) * t36) * mrSges(6,3) + m(6) * (t10 * t3 + t20 * t19 + t9 * t2) + m(5) * (t17 * t12 + t18 * t13 + t49 * t41) + (t48 - pkin(2)) * t42 + (t49 + t41) * t15 + (t19 + t20) * t14; -0.2e1 * pkin(2) * t42 + t20 * t91 + t49 * t90 + (t10 * t36 - t9 * t37) * t74 + (-t17 * t37 + t18 * t36) * t75 + pkin(7) * t68 + m(6) * (t10 ^ 2 + t20 ^ 2 + t9 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2 + t49 ^ 2) + m(4) * (t77 * pkin(7) ^ 2 + pkin(2) ^ 2) + t67; t65 + m(6) * (t47 * t2 + t3 * t87) + t69 * t46 + (t12 * t59 + t13 * t56) * t89 + t94; t65 + (t17 * t59 + t18 * t56) * t89 + m(6) * (t10 * t87 + t47 * t9) + t69 * pkin(7) + t93; t47 * t92 + Ifges(4,3) + m(6) * (t47 ^ 2 + t52) + m(5) * (t59 ^ 2 * t63 + t52) + 0.2e1 * t95 + t78; (m(6) * t2 - t80) * pkin(4) + t70 + t94; (m(6) * t9 - t80) * pkin(4) + t70 + t93; t47 * t88 + (pkin(4) + t47) * mrSges(6,1) + t95 + t78; (t92 + t88) * pkin(4) + t78; m(6) * t19 + t14; m(6) * t20 + t14; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
