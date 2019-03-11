% Calculate joint inertia matrix for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:59
% EndTime: 2019-03-08 18:22:59
% DurationCPUTime: 0.07s
% Computational Cost: add. (38->28), mult. (64->31), div. (0->0), fcn. (15->2), ass. (0->8)
t8 = 2 * mrSges(5,3);
t7 = Ifges(5,2) + Ifges(4,3);
t3 = sin(qJ(3));
t4 = cos(qJ(3));
t6 = (mrSges(4,1) * t4 - mrSges(4,2) * t3) * pkin(2);
t2 = -pkin(2) * t4 - pkin(3);
t1 = pkin(2) * t3 + qJ(4);
t5 = [m(2) + m(3) + m(4) + m(5); 0; -0.2e1 * t2 * mrSges(5,1) + t1 * t8 + Ifges(3,3) + 0.2e1 * t6 + m(5) * (t1 ^ 2 + t2 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) * pkin(2) ^ 2 + t7; 0; m(5) * (-pkin(3) * t2 + qJ(4) * t1) + t6 + (qJ(4) + t1) * mrSges(5,3) + (pkin(3) - t2) * mrSges(5,1) + t7; 0.2e1 * pkin(3) * mrSges(5,1) + qJ(4) * t8 + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t7; 0; m(5) * t2 - mrSges(5,1); -m(5) * pkin(3) - mrSges(5,1); m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1) t5(2) t5(4) t5(7); t5(2) t5(3) t5(5) t5(8); t5(4) t5(5) t5(6) t5(9); t5(7) t5(8) t5(9) t5(10);];
Mq  = res;
