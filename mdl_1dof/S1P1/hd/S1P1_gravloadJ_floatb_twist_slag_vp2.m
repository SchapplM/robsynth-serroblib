% Calculate Gravitation load on the joints for
% S1P1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m [2x1]
%   mass of all robot links (including the base)
% mrSges [2x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S1P1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(1,1),zeros(2,1),zeros(2,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S1P1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [2 1]), ...
  'S1P1_gravloadJ_floatb_twist_slag_vp2: m has to be [2x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [2,3]), ...
  'S1P1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [2x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:52
% EndTime: 2021-01-14 12:21:53
% DurationCPUTime: 0.13s
% Computational Cost: add. (1->1), mult. (1->1), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [-g(3) * m(2)];
taug = t1(:);
