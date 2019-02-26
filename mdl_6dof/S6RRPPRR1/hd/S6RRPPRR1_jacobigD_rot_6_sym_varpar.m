% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:07
% EndTime: 2019-02-26 21:28:07
% DurationCPUTime: 0.03s
% Computational Cost: add. (31->8), mult. (48->12), div. (0->0), fcn. (48->6), ass. (0->13)
t155 = -qJD(2) + qJD(5);
t147 = sin(qJ(1));
t154 = qJD(1) * t147;
t149 = cos(qJ(1));
t153 = qJD(1) * t149;
t145 = qJ(2) + pkin(10);
t143 = sin(t145);
t144 = cos(t145);
t146 = sin(qJ(5));
t148 = cos(qJ(5));
t152 = t143 * t148 - t144 * t146;
t150 = t155 * (t143 * t146 + t144 * t148);
t1 = [0, t153, 0, 0, -t153, t150 * t149 + t152 * t154; 0, t154, 0, 0, -t154, t150 * t147 - t152 * t153; 0, 0, 0, 0, 0, t155 * t152;];
JgD_rot  = t1;
