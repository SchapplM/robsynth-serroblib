% Calculate vector of inverse dynamics joint torques for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% qJDD [7x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% mrSges [8x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [8x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [7x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 19:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S7RRRRRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1),zeros(8,1),zeros(8,3),zeros(8,6)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_slag_vp2: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_slag_vp2: qJD has to be [7x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [7 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_slag_vp2: m has to be [8x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [8,3]), ...
  'S7RRRRRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [8x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [8 6]), ...
  'S7RRRRRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [8x6] (double)'); 

%% Function calls and calculation
g = S7RRRRRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges);
c = S7RRRRRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges);
M = S7RRRRRRR1_inertiaJ_slag_vp2(qJ, ...
    pkin, m, mrSges, Ifges);
tau = g+c+M*qJDD;

tau  = tau;
