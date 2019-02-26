% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR2_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:15
% EndTime: 2019-02-26 20:19:15
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->12), mult. (60->30), div. (0->0), fcn. (62->8), ass. (0->21)
t171 = qJ(3) + qJ(4);
t168 = sin(t171);
t173 = sin(pkin(6));
t184 = t173 * t168;
t175 = cos(pkin(6));
t176 = sin(qJ(2));
t183 = t175 * t176;
t177 = cos(qJ(2));
t182 = t175 * t177;
t181 = qJD(2) * t168;
t180 = qJD(2) * t173;
t172 = sin(pkin(12));
t174 = cos(pkin(12));
t179 = t172 * t177 + t174 * t183;
t178 = -t172 * t183 + t174 * t177;
t170 = qJD(3) + qJD(4);
t169 = cos(t171);
t167 = t176 * t180;
t166 = t178 * qJD(2);
t165 = t179 * qJD(2);
t1 = [0, 0, t166, t166 (t178 * t169 + t172 * t184) * t170 + (-t172 * t182 - t174 * t176) * t181, 0; 0, 0, t165, t165 (t179 * t169 - t174 * t184) * t170 + (-t172 * t176 + t174 * t182) * t181, 0; 0, 0, t167, t167, t173 * t176 * t170 * t169 + (t170 * t175 + t177 * t180) * t168, 0;];
JgD_rot  = t1;
