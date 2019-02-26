% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRP4_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobigD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:59
% EndTime: 2019-02-26 20:16:59
% DurationCPUTime: 0.04s
% Computational Cost: add. (22->10), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
t154 = sin(pkin(6));
t157 = sin(qJ(3));
t167 = t154 * t157;
t158 = sin(qJ(2));
t166 = t154 * t158;
t156 = cos(pkin(6));
t165 = t156 * t158;
t160 = cos(qJ(2));
t164 = t156 * t160;
t163 = qJD(2) * t157;
t153 = sin(pkin(11));
t155 = cos(pkin(11));
t162 = t153 * t160 + t155 * t165;
t161 = -t153 * t165 + t155 * t160;
t159 = cos(qJ(3));
t152 = t154 * t160 * t163 + (t156 * t157 + t159 * t166) * qJD(3);
t151 = (t153 * t167 + t161 * t159) * qJD(3) + (-t153 * t164 - t155 * t158) * t163;
t150 = (-t155 * t167 + t162 * t159) * qJD(3) + (-t153 * t158 + t155 * t164) * t163;
t1 = [0, 0, t161 * qJD(2), t151, t151, 0; 0, 0, t162 * qJD(2), t150, t150, 0; 0, 0, qJD(2) * t166, t152, t152, 0;];
JgD_rot  = t1;
