% Calculate vector of cutting forces with Newton-Euler
% S1P1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
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
% Ifges [2x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x2]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S1P1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(1,1),zeros(2,1),zeros(2,3),zeros(2,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_invdynf_fixb_snew_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1P1_invdynf_fixb_snew_vp2: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'S1P1_invdynf_fixb_snew_vp2: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S1P1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_invdynf_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [2 1]), ...
  'S1P1_invdynf_fixb_snew_vp2: m has to be [2x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [2,3]), ...
  'S1P1_invdynf_fixb_snew_vp2: mrSges has to be [2x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [2 6]), ...
  'S1P1_invdynf_fixb_snew_vp2: Ifges has to be [2x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:53
% EndTime: 2021-01-14 12:21:53
% DurationCPUTime: 0.07s
% Computational Cost: add. (7->5), mult. (9->6), div. (0->0), fcn. (0->0), ass. (0->3)
t2 = -m(1) - m(2);
t1 = m(2) * (-g(3) + qJDD(1));
t3 = [t2 * g(1), -m(2) * g(1); t2 * g(2), -m(2) * g(2); -m(1) * g(3) + t1, t1;];
f_new = t3;
