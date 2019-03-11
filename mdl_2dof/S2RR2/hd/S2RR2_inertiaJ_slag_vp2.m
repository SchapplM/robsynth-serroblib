% Calculate joint inertia matrix for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S2RR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_inertiaJ_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_inertiaJ_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR2_inertiaJ_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR2_inertiaJ_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR2_inertiaJ_slag_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:52
% EndTime: 2019-03-08 18:00:53
% DurationCPUTime: 0.06s
% Computational Cost: add. (11->9), mult. (30->16), div. (0->0), fcn. (14->2), ass. (0->4)
t4 = cos(qJ(2));
t2 = t4 ^ 2;
t3 = sin(qJ(2));
t1 = [Ifges(3,2) * t2 + Ifges(2,3) + (Ifges(3,1) * t3 + 0.2e1 * Ifges(3,4) * t4) * t3 + (m(3) * pkin(1) + 2 * mrSges(3,3)) * (t3 ^ 2 + t2) * pkin(1); Ifges(3,5) * t3 + Ifges(3,6) * t4 + (-mrSges(3,1) * t3 - mrSges(3,2) * t4) * pkin(1); Ifges(3,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t1(1) t1(2); t1(2) t1(3);];
Mq  = res;
