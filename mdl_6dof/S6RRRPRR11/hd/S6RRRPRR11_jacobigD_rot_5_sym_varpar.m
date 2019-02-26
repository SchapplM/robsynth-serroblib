% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR11_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobigD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:50
% EndTime: 2019-02-26 22:21:50
% DurationCPUTime: 0.09s
% Computational Cost: add. (13->9), mult. (48->16), div. (0->0), fcn. (48->6), ass. (0->15)
t147 = sin(qJ(2));
t148 = sin(qJ(1));
t156 = t147 * t148;
t150 = cos(qJ(1));
t155 = t147 * t150;
t149 = cos(qJ(2));
t154 = t148 * t149;
t153 = t149 * t150;
t145 = sin(pkin(6));
t152 = qJD(1) * t145;
t151 = t145 * qJD(2) * t147;
t146 = cos(pkin(6));
t144 = (t146 * t155 + t154) * qJD(2) + (t146 * t154 + t155) * qJD(1);
t143 = (t146 * t156 - t153) * qJD(2) + (-t146 * t153 + t156) * qJD(1);
t1 = [0, t150 * t152, -t143, 0, t143, 0; 0, t148 * t152, t144, 0, -t144, 0; 0, 0, t151, 0, -t151, 0;];
JgD_rot  = t1;
