% Calculate joint inertia matrix for
% S4PPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:05
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PPPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR5_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:01
% EndTime: 2018-11-14 14:05:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (20->15), mult. (47->24), div. (0->0), fcn. (35->4), ass. (0->7)
t7 = cos(qJ(4));
t6 = sin(qJ(4));
t5 = cos(pkin(5));
t4 = sin(pkin(5));
t2 = t4 * t7 - t5 * t6;
t1 = -t4 * t6 - t5 * t7;
t3 = [m(2) + m(5) * (t1 ^ 2 + t2 ^ 2) + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * (t4 ^ 2 + t5 ^ 2); 0; m(3) + m(4) + m(5); -m(4) * t5 + m(5) * (t1 * t7 + t2 * t6); 0; m(4) + m(5) * (t6 ^ 2 + t7 ^ 2); mrSges(5,1) * t1 - mrSges(5,2) * t2; 0; mrSges(5,1) * t7 - mrSges(5,2) * t6; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1) t3(2) t3(4) t3(7); t3(2) t3(3) t3(5) t3(8); t3(4) t3(5) t3(6) t3(9); t3(7) t3(8) t3(9) t3(10);];
Mq  = res;
