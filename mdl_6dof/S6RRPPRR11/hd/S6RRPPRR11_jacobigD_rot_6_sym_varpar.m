% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR11_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:20
% EndTime: 2019-02-26 21:34:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
t166 = sin(pkin(6));
t169 = sin(qJ(1));
t184 = t166 * t169;
t171 = cos(qJ(1));
t183 = t166 * t171;
t168 = sin(qJ(2));
t182 = t169 * t168;
t170 = cos(qJ(2));
t181 = t169 * t170;
t180 = t170 * t171;
t179 = t171 * t168;
t178 = qJD(1) * t166;
t165 = pkin(11) + qJ(5);
t164 = cos(t165);
t177 = qJD(2) * t164;
t176 = qJD(2) * t166;
t167 = cos(pkin(6));
t175 = t167 * t180 - t182;
t174 = t167 * t181 + t179;
t173 = t167 * t179 + t181;
t172 = -t167 * t182 + t180;
t163 = sin(t165);
t1 = [0, t171 * t178, 0, 0, -t173 * qJD(1) - t174 * qJD(2) (t174 * t163 + t164 * t184) * qJD(5) - t172 * t177 + (t163 * t183 - t175 * t164) * qJD(1); 0, t169 * t178, 0, 0, t172 * qJD(1) + t175 * qJD(2) (-t175 * t163 - t164 * t183) * qJD(5) - t173 * t177 + (t163 * t184 - t174 * t164) * qJD(1); 0, 0, 0, 0, t170 * t176, -t168 * t164 * t176 + (-t163 * t166 * t170 + t164 * t167) * qJD(5);];
JgD_rot  = t1;
