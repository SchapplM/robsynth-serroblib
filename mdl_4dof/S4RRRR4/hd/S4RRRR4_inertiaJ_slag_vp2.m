% Calculate joint inertia matrix for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:41
% EndTime: 2019-12-31 17:25:42
% DurationCPUTime: 0.32s
% Computational Cost: add. (421->109), mult. (843->169), div. (0->0), fcn. (764->6), ass. (0->49)
t49 = cos(qJ(2));
t73 = -pkin(6) - pkin(5);
t32 = t73 * t49;
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t46 = sin(qJ(2));
t59 = t73 * t46;
t16 = -t45 * t32 - t48 * t59;
t27 = t45 * t49 + t48 * t46;
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t56 = mrSges(5,1) * t44 + mrSges(5,2) * t47;
t9 = t56 * t27;
t78 = m(5) * t16 + t9;
t26 = t45 * t46 - t48 * t49;
t67 = t44 * mrSges(5,3);
t11 = -t26 * mrSges(5,2) - t27 * t67;
t68 = t27 * t47;
t12 = t26 * mrSges(5,1) - mrSges(5,3) * t68;
t37 = -t49 * pkin(2) - pkin(1);
t10 = t26 * pkin(3) - t27 * pkin(7) + t37;
t18 = -t48 * t32 + t45 * t59;
t2 = t47 * t10 - t44 * t18;
t3 = t44 * t10 + t47 * t18;
t72 = t3 * t47;
t77 = m(5) * (-t2 * t44 + t72) + t47 * t11 - t44 * t12;
t76 = t16 ^ 2;
t75 = 0.2e1 * t16;
t74 = 0.2e1 * t37;
t71 = Ifges(5,4) * t44;
t70 = Ifges(5,4) * t47;
t69 = t27 * t44;
t64 = Ifges(5,5) * t68 + Ifges(5,3) * t26;
t63 = Ifges(5,5) * t44 + Ifges(5,6) * t47;
t62 = t44 ^ 2 + t47 ^ 2;
t61 = t46 ^ 2 + t49 ^ 2;
t30 = Ifges(5,2) * t47 + t71;
t31 = Ifges(5,1) * t44 + t70;
t60 = t47 * t30 + t44 * t31 + Ifges(4,3);
t35 = t45 * pkin(2) + pkin(7);
t58 = t62 * t35;
t55 = 0.2e1 * t62 * mrSges(5,3);
t54 = (t48 * mrSges(4,1) - t45 * mrSges(4,2)) * pkin(2);
t29 = -t47 * mrSges(5,1) + t44 * mrSges(5,2);
t6 = Ifges(5,6) * t26 + (-Ifges(5,2) * t44 + t70) * t27;
t7 = Ifges(5,5) * t26 + (Ifges(5,1) * t47 - t71) * t27;
t53 = -t18 * mrSges(4,2) + mrSges(5,3) * t72 - t2 * t67 + t44 * t7 / 0.2e1 - t30 * t69 / 0.2e1 + t31 * t68 / 0.2e1 + Ifges(4,5) * t27 + t47 * t6 / 0.2e1 + (t63 / 0.2e1 - Ifges(4,6)) * t26 + (-mrSges(4,1) + t29) * t16;
t36 = -t48 * pkin(2) - pkin(3);
t1 = [Ifges(2,3) + t9 * t75 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 - 0.2e1 * pkin(1) * (-t49 * mrSges(3,1) + t46 * mrSges(3,2)) + t46 * (Ifges(3,1) * t46 + Ifges(3,4) * t49) + t49 * (Ifges(3,4) * t46 + Ifges(3,2) * t49) + 0.2e1 * t61 * pkin(5) * mrSges(3,3) + (mrSges(4,1) * t74 - 0.2e1 * t18 * mrSges(4,3) + Ifges(4,2) * t26 + t64) * t26 + (mrSges(4,2) * t74 + mrSges(4,3) * t75 + Ifges(4,1) * t27 - t44 * t6 + t47 * t7 + (-Ifges(5,6) * t44 - (2 * Ifges(4,4))) * t26) * t27 + m(5) * (t2 ^ 2 + t3 ^ 2 + t76) + m(4) * (t18 ^ 2 + t37 ^ 2 + t76) + m(3) * (t61 * pkin(5) ^ 2 + pkin(1) ^ 2); t53 + (m(4) * (-t16 * t48 + t18 * t45) + (-t45 * t26 - t48 * t27) * mrSges(4,3)) * pkin(2) + (-t46 * mrSges(3,1) - t49 * mrSges(3,2)) * pkin(5) + Ifges(3,5) * t46 + Ifges(3,6) * t49 + t78 * t36 + t77 * t35; 0.2e1 * t36 * t29 + Ifges(3,3) + 0.2e1 * t54 + t35 * t55 + m(5) * (t62 * t35 ^ 2 + t36 ^ 2) + m(4) * (t45 ^ 2 + t48 ^ 2) * pkin(2) ^ 2 + t60; -t78 * pkin(3) + t77 * pkin(7) + t53; m(5) * (-pkin(3) * t36 + pkin(7) * t58) + (-pkin(3) + t36) * t29 + t54 + (t62 * pkin(7) + t58) * mrSges(5,3) + t60; -0.2e1 * pkin(3) * t29 + m(5) * (t62 * pkin(7) ^ 2 + pkin(3) ^ 2) + pkin(7) * t55 + t60; t2 * mrSges(5,1) - t3 * mrSges(5,2) - Ifges(5,6) * t69 + t64; -t56 * t35 + t63; -t56 * pkin(7) + t63; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
