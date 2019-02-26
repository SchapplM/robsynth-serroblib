% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRRR3_jacobigD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobigD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobigD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobigD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:05
% EndTime: 2019-02-26 19:44:05
% DurationCPUTime: 0.04s
% Computational Cost: add. (13->13), mult. (43->29), div. (0->0), fcn. (47->11), ass. (0->15)
t159 = sin(pkin(13));
t166 = cos(pkin(6));
t172 = t159 * t166;
t161 = sin(pkin(7));
t162 = sin(pkin(6));
t171 = t162 * t161;
t164 = cos(pkin(13));
t170 = t164 * t166;
t169 = qJD(3) * sin(pkin(8));
t168 = cos(qJ(3));
t167 = sin(qJ(3));
t165 = cos(pkin(7));
t163 = cos(pkin(14));
t158 = sin(pkin(14));
t1 = [0, 0, 0 -(-(-t158 * t172 + t164 * t163) * t168 + (-(-t164 * t158 - t163 * t172) * t165 - t159 * t171) * t167) * t169, 0, 0; 0, 0, 0 -(-(t158 * t170 + t159 * t163) * t168 + (-(-t159 * t158 + t163 * t170) * t165 + t164 * t171) * t167) * t169, 0, 0; 0, 0, 0 -(-t161 * t166 * t167 + (-t163 * t165 * t167 - t158 * t168) * t162) * t169, 0, 0;];
JgD_rot  = t1;
