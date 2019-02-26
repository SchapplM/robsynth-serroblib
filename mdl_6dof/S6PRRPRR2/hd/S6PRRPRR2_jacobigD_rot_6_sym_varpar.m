% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:42
% EndTime: 2019-02-26 20:04:42
% DurationCPUTime: 0.04s
% Computational Cost: add. (40->11), mult. (84->30), div. (0->0), fcn. (88->8), ass. (0->20)
t167 = qJ(3) + pkin(12);
t165 = sin(t167);
t169 = sin(pkin(6));
t180 = t169 * t165;
t171 = cos(pkin(6));
t172 = sin(qJ(2));
t179 = t171 * t172;
t173 = cos(qJ(2));
t178 = t171 * t173;
t177 = qJD(2) * t165;
t176 = qJD(2) * t169;
t168 = sin(pkin(11));
t170 = cos(pkin(11));
t175 = t168 * t173 + t170 * t179;
t174 = -t168 * t179 + t170 * t173;
t166 = cos(t167);
t164 = t173 * t165 * t176 + (t166 * t169 * t172 + t165 * t171) * qJD(3);
t163 = (t174 * t166 + t168 * t180) * qJD(3) + (-t168 * t178 - t170 * t172) * t177;
t162 = (t175 * t166 - t170 * t180) * qJD(3) + (-t168 * t172 + t170 * t178) * t177;
t1 = [0, 0, t174 * qJD(2), 0, t163, t163; 0, 0, t175 * qJD(2), 0, t162, t162; 0, 0, t172 * t176, 0, t164, t164;];
JgD_rot  = t1;
