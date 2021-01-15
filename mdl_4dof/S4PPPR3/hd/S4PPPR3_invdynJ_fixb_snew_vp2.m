% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta3]';
% m [5x1]
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PPPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:23
% EndTime: 2021-01-15 14:39:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (119->23), mult. (124->31), div. (0->0), fcn. (80->4), ass. (0->16)
t107 = qJD(4) ^ 2;
t106 = cos(qJ(4));
t105 = sin(qJ(4));
t104 = cos(pkin(5));
t103 = sin(pkin(5));
t102 = -g(1) + qJDD(2);
t101 = -g(2) + qJDD(1);
t99 = t104 * t101 + t103 * t102;
t98 = -t103 * t101 + t104 * t102;
t97 = t105 * t98 + t106 * t99;
t96 = -t105 * t99 + t106 * t98;
t95 = m(5) * t97 - t107 * mrSges(5,1) - qJDD(4) * mrSges(5,2);
t94 = m(5) * t96 + qJDD(4) * mrSges(5,1) - t107 * mrSges(5,2);
t93 = m(4) * t99 - t105 * t94 + t106 * t95;
t92 = m(4) * t98 + t105 * t95 + t106 * t94;
t1 = [-t103 * t92 + t104 * t93 + (m(2) + m(3)) * t101; m(3) * t102 + t103 * t93 + t104 * t92; (m(4) + m(5)) * (g(3) + qJDD(3)); mrSges(5,1) * t96 - mrSges(5,2) * t97 + Ifges(5,3) * qJDD(4);];
tauJ = t1;
