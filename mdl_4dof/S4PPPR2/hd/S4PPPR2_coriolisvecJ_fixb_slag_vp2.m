% Calculate vector of centrifugal and coriolis load on the joints for
% S4PPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% tauc [4x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:56
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PPPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:55:24
% EndTime: 2018-11-14 13:55:24
% DurationCPUTime: 0.09s
% Computational Cost: add. (53->15), mult. (140->32), div. (0->0), fcn. (94->4), ass. (0->15)
t8 = sin(pkin(5));
t16 = qJD(1) * t8;
t10 = sin(qJ(4));
t11 = cos(qJ(4));
t9 = cos(pkin(5));
t7 = -qJD(1) * t9 + qJD(3);
t3 = -t10 * t16 + t11 * t7;
t1 = qJD(4) * t3;
t4 = t10 * t7 + t11 * t16;
t15 = qJD(4) * t4;
t14 = t10 * t9 - t11 * t8;
t13 = -t10 * t8 - t11 * t9;
t6 = t13 * qJD(4);
t5 = t14 * qJD(4);
t2 = [m(5) * (-t1 * t14 - t13 * t15 + t3 * t5 + t4 * t6) + (t5 * mrSges(5,1) - t6 * mrSges(5,2)) * qJD(4); 0; m(5) * (t1 * t10 - t11 * t15 + (-t10 * t3 + t11 * t4) * qJD(4)) + (-t10 * mrSges(5,1) - t11 * mrSges(5,2)) * qJD(4) ^ 2; 0;];
tauc  = t2(:);
