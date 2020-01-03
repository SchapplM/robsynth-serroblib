% Calculate joint inertia matrix for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:03
% DurationCPUTime: 0.14s
% Computational Cost: add. (87->48), mult. (160->54), div. (0->0), fcn. (61->2), ass. (0->14)
t10 = cos(qJ(3));
t9 = sin(qJ(3));
t16 = t10 ^ 2 + t9 ^ 2;
t22 = (-mrSges(4,3) - mrSges(5,2)) * t16;
t20 = -m(5) * pkin(3) - mrSges(5,1);
t19 = 2 * mrSges(5,1);
t18 = 2 * qJ(2);
t11 = -pkin(1) - pkin(5);
t17 = t16 * t11 ^ 2;
t14 = 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t16;
t13 = (m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t9 + (mrSges(4,1) - t20) * t10;
t12 = qJ(2) ^ 2;
t2 = t9 * pkin(3) - t10 * qJ(4) + qJ(2);
t1 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t18) + Ifges(3,1) + Ifges(2,3) + (mrSges(4,1) * t18 + t2 * t19 + (Ifges(5,3) + Ifges(4,2)) * t9) * t9 + (mrSges(4,2) * t18 - 0.2e1 * t2 * mrSges(5,3) + (Ifges(5,1) + Ifges(4,1)) * t10 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t9) * t10 + m(5) * (t2 ^ 2 + t17) + m(4) * (t12 + t17) + (m(3) * (pkin(1) ^ 2 + t12)) + 0.2e1 * t11 * t22; -(m(3) * pkin(1)) + t11 * t14 + mrSges(3,2) + t22; m(3) + t14; (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t9 + (-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t10 + t13 * t11; t13; Ifges(5,2) + Ifges(4,3) + pkin(3) * t19 + 0.2e1 * qJ(4) * mrSges(5,3) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2); (-m(5) * t11 + mrSges(5,2)) * t10; -m(5) * t10; t20; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
