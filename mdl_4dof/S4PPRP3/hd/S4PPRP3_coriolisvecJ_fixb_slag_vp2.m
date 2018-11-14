% Calculate vector of centrifugal and coriolis load on the joints for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
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
% tauc [4x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:58
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PPRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:58:14
% EndTime: 2018-11-14 13:58:15
% DurationCPUTime: 0.14s
% Computational Cost: add. (70->22), mult. (193->34), div. (0->0), fcn. (90->2), ass. (0->15)
t12 = sin(qJ(3));
t13 = cos(qJ(3));
t9 = -t12 * qJD(1) + t13 * qJD(2);
t7 = qJD(3) * t9;
t10 = t13 * qJD(1) + t12 * qJD(2);
t8 = qJD(3) * t10;
t21 = t10 * t12;
t20 = mrSges(4,1) + mrSges(5,1);
t19 = mrSges(4,2) + mrSges(5,2);
t18 = qJD(3) * t12;
t16 = t7 * t12;
t15 = t8 * t12 + t7 * t13;
t14 = qJD(3) ^ 2;
t6 = qJD(3) * pkin(3) + t9;
t1 = [m(4) * ((-t13 * t9 - t21) * qJD(3) + t15) + m(5) * ((-t13 * t6 - t21) * qJD(3) + t15) + (t19 * t12 - t20 * t13) * t14; m(4) * (-t9 * t18 + t16) + m(5) * (-t6 * t18 + t16) + (-t20 * t12 - t19 * t13) * t14; -(m(5) * pkin(3) + t20) * t8 + (-m(5) * (-t6 + t9) + t20 * qJD(3)) * t10; 0;];
tauc  = t1(:);
