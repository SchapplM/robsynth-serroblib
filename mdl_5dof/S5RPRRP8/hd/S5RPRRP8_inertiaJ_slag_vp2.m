% Calculate joint inertia matrix for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:08
% EndTime: 2019-12-31 18:47:09
% DurationCPUTime: 0.34s
% Computational Cost: add. (318->94), mult. (552->110), div. (0->0), fcn. (310->4), ass. (0->33)
t37 = cos(qJ(4));
t33 = t37 ^ 2;
t35 = sin(qJ(4));
t49 = t35 ^ 2 + t33;
t15 = t37 * mrSges(6,1) + t35 * mrSges(6,3);
t16 = -t37 * mrSges(5,1) + t35 * mrSges(5,2);
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t64 = (mrSges(6,2) + mrSges(5,3)) * t49;
t71 = -(mrSges(4,2) - t64) * t36 + (mrSges(4,1) + t15 - t16) * t38;
t68 = t49 * t36;
t65 = -m(6) * pkin(4) - mrSges(6,1);
t63 = -(pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t35 + (qJ(5) * mrSges(6,2) + Ifges(5,6) - Ifges(6,6)) * t37;
t42 = -(-Ifges(5,2) - Ifges(6,3)) * t33 + Ifges(4,3) + (0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t37 + (Ifges(5,1) + Ifges(6,1)) * t35) * t35;
t62 = -0.2e1 * t16;
t39 = -pkin(1) - pkin(2);
t13 = t38 * qJ(2) + t36 * t39;
t11 = -pkin(7) + t13;
t61 = t68 * t11;
t60 = t49 * pkin(7) * t11;
t59 = t49 * t11 ^ 2;
t12 = -t36 * qJ(2) + t38 * t39;
t57 = t12 * mrSges(4,1);
t56 = t13 * mrSges(4,2);
t51 = t68 * pkin(7);
t50 = t49 * pkin(7) ^ 2;
t14 = -t37 * pkin(4) - t35 * qJ(5) - pkin(3);
t41 = (m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t37 + (-mrSges(5,1) + t65) * t35;
t34 = t38 ^ 2;
t32 = t36 ^ 2;
t10 = pkin(3) - t12;
t1 = -t12 - t14;
t2 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t57 + 0.2e1 * t56 + 0.2e1 * qJ(2) * mrSges(3,3) + 0.2e1 * t1 * t15 + t10 * t62 + Ifges(3,2) + Ifges(2,3) + m(5) * (t10 ^ 2 + t59) + m(6) * (t1 ^ 2 + t59) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) + t42 - 0.2e1 * t11 * t64; -m(3) * pkin(1) - mrSges(3,1) + m(5) * (-t38 * t10 + t61) + m(6) * (-t38 * t1 + t61) + m(4) * (t38 * t12 + t36 * t13) - t71; m(3) + m(4) * (t32 + t34) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t49 * t32 + t34); t57 - t56 + (pkin(3) + t10) * t16 + (t14 - t1) * t15 + m(5) * (-pkin(3) * t10 + t60) + m(6) * (t14 * t1 + t60) + (-pkin(7) + t11) * t64 - t42; m(5) * (pkin(3) * t38 + t51) + m(6) * (-t14 * t38 + t51) + t71; pkin(3) * t62 - 0.2e1 * t14 * t15 + m(5) * (pkin(3) ^ 2 + t50) + m(6) * (t14 ^ 2 + t50) + t42 + 0.2e1 * pkin(7) * t64; t41 * t11 - t63; t41 * t36; t41 * pkin(7) + t63; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); (m(6) * t11 - mrSges(6,2)) * t35; m(6) * t35 * t36; (m(6) * pkin(7) + mrSges(6,2)) * t35; t65; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
