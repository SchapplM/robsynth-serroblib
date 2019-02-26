% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR9_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:42:23
% EndTime: 2019-02-26 21:42:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
t169 = sin(pkin(6));
t172 = sin(qJ(1));
t187 = t169 * t172;
t174 = cos(qJ(1));
t186 = t169 * t174;
t171 = sin(qJ(2));
t185 = t171 * t172;
t184 = t171 * t174;
t173 = cos(qJ(2));
t183 = t172 * t173;
t182 = t174 * t173;
t181 = qJD(1) * t169;
t168 = pkin(11) + qJ(4);
t166 = sin(t168);
t180 = qJD(2) * t166;
t179 = qJD(2) * t169;
t170 = cos(pkin(6));
t178 = t170 * t182 - t185;
t177 = t170 * t183 + t184;
t176 = t170 * t184 + t183;
t175 = -t170 * t185 + t182;
t167 = cos(t168);
t1 = [0, t174 * t181, 0, t178 * qJD(1) + t175 * qJD(2), 0 (t166 * t187 + t175 * t167) * qJD(4) - t177 * t180 + (-t176 * t166 - t167 * t186) * qJD(1); 0, t172 * t181, 0, t177 * qJD(1) + t176 * qJD(2), 0 (-t166 * t186 + t176 * t167) * qJD(4) + t178 * t180 + (t175 * t166 - t167 * t187) * qJD(1); 0, 0, 0, t171 * t179, 0, t173 * t166 * t179 + (t167 * t169 * t171 + t166 * t170) * qJD(4);];
JgD_rot  = t1;
