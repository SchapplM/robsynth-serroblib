% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRP3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:23
% EndTime: 2019-02-26 20:02:23
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t162 = sin(pkin(6));
t165 = sin(qJ(3));
t175 = t162 * t165;
t166 = sin(qJ(2));
t174 = t162 * t166;
t164 = cos(pkin(6));
t173 = t164 * t166;
t168 = cos(qJ(2));
t172 = t164 * t168;
t171 = qJD(2) * t165;
t161 = sin(pkin(10));
t163 = cos(pkin(10));
t170 = t161 * t168 + t163 * t173;
t169 = -t161 * t173 + t163 * t168;
t167 = cos(qJ(3));
t1 = [0, 0, t169 * qJD(2), 0 (t161 * t175 + t169 * t167) * qJD(3) + (-t161 * t172 - t163 * t166) * t171, 0; 0, 0, t170 * qJD(2), 0 (-t163 * t175 + t170 * t167) * qJD(3) + (-t161 * t166 + t163 * t172) * t171, 0; 0, 0, qJD(2) * t174, 0, t162 * t168 * t171 + (t164 * t165 + t167 * t174) * qJD(3), 0;];
JgD_rot  = t1;
