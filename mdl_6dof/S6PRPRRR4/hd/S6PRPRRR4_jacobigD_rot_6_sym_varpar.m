% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:55:23
% EndTime: 2019-02-26 19:55:23
% DurationCPUTime: 0.04s
% Computational Cost: add. (40->11), mult. (84->30), div. (0->0), fcn. (88->8), ass. (0->20)
t164 = pkin(12) + qJ(4);
t162 = sin(t164);
t166 = sin(pkin(6));
t177 = t166 * t162;
t168 = cos(pkin(6));
t169 = sin(qJ(2));
t176 = t168 * t169;
t170 = cos(qJ(2));
t175 = t168 * t170;
t174 = qJD(2) * t162;
t173 = qJD(2) * t166;
t165 = sin(pkin(11));
t167 = cos(pkin(11));
t172 = t165 * t170 + t167 * t176;
t171 = -t165 * t176 + t167 * t170;
t163 = cos(t164);
t161 = t170 * t162 * t173 + (t163 * t166 * t169 + t162 * t168) * qJD(4);
t160 = (t171 * t163 + t165 * t177) * qJD(4) + (-t165 * t175 - t167 * t169) * t174;
t159 = (t172 * t163 - t167 * t177) * qJD(4) + (-t165 * t169 + t167 * t175) * t174;
t1 = [0, 0, 0, t171 * qJD(2), t160, t160; 0, 0, 0, t172 * qJD(2), t159, t159; 0, 0, 0, t169 * t173, t161, t161;];
JgD_rot  = t1;
