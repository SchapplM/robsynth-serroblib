% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR3
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
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:49
% EndTime: 2019-02-26 19:54:49
% DurationCPUTime: 0.04s
% Computational Cost: add. (38->12), mult. (60->30), div. (0->0), fcn. (62->8), ass. (0->21)
t171 = pkin(12) + qJ(4) + qJ(5);
t169 = sin(t171);
t174 = sin(pkin(6));
t185 = t174 * t169;
t176 = cos(pkin(6));
t177 = sin(qJ(2));
t184 = t176 * t177;
t178 = cos(qJ(2));
t183 = t176 * t178;
t182 = qJD(2) * t169;
t181 = qJD(2) * t174;
t173 = sin(pkin(11));
t175 = cos(pkin(11));
t180 = t173 * t178 + t175 * t184;
t179 = -t173 * t184 + t175 * t178;
t172 = qJD(4) + qJD(5);
t170 = cos(t171);
t168 = t177 * t181;
t167 = t179 * qJD(2);
t166 = t180 * qJD(2);
t1 = [0, 0, 0, t167, t167 (t179 * t170 + t173 * t185) * t172 + (-t173 * t183 - t175 * t177) * t182; 0, 0, 0, t166, t166 (t180 * t170 - t175 * t185) * t172 + (-t173 * t177 + t175 * t183) * t182; 0, 0, 0, t168, t168, t174 * t177 * t172 * t170 + (t172 * t176 + t178 * t181) * t169;];
JgD_rot  = t1;
